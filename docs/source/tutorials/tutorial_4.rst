Tutorial #4: File Reading and Antiferromagnets
==============================================

**Aims**:

* Set up a structure using input files rather than by hand
* Work with antiferromagnetic structures
* Extract raw bandstructure data

In this example, we'll use an iron-gadolinium toy model:

.. raw:: html

    <div id="viewer" style="width: 100%; max-width: 400px; height: 400px; border: 1px solid #ccc;"></div>

    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    window.onload = () => {
      const viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });

      const a = 10.0;
      const numRepeats = 1;
      let atomLines = [];

      for (let i = 0; i <= numRepeats; i++) {
        for (let j = 0; j <= numRepeats; j++) {
          for (let k = 0; k <= numRepeats; k++) {
            const x = i * a, y = j * a, z = k * a;
            atomLines.push(`Fe ${x.toFixed(3)} ${y.toFixed(3)} ${z.toFixed(3)}`);
          }
        }
      }

      for (let i = 0; i < numRepeats; i++) {
        for (let j = 0; j < numRepeats; j++) {
          for (let k = 0; k < numRepeats; k++) {
            const x = (i + 0.5) * a;
            const y = (j + 0.5) * a;
            const z = (k + 0.5) * a;
            atomLines.push(`Gd ${x.toFixed(3)} ${y.toFixed(3)} ${z.toFixed(3)}`);
          }
        }
      }

    const xyz = `${atomLines.length}
    Fe corners + Gd centers
    ${atomLines.join("\n")}`;
      const model = viewer.addModel(xyz, "xyz");

      model.setStyle({ elem: "Fe" }, { sphere: { scale: 0.3, color: "gold" }, stick: {} });
      model.setStyle({ elem: "Gd" }, { sphere: { scale: 0.5, color: "purple" }, stick: {} });

      const drawUnitCell = (x0, y0, z0, a) => {
        const corners = [
          [x0, y0, z0], [x0 + a, y0, z0], [x0, y0 + a, z0], [x0 + a, y0 + a, z0],
          [x0, y0, z0 + a], [x0 + a, y0, z0 + a], [x0, y0 + a, z0 + a], [x0 + a, y0 + a, z0 + a]
        ];

        const edges = [
          [0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6],
          [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]
        ];

        for (const [i, j] of edges) {
          const p1 = corners[i];
          const p2 = corners[j];
          viewer.addLine({
            start: { x: p1[0], y: p1[1], z: p1[2] },
            end: { x: p2[0], y: p2[1], z: p2[2] },
            color: "black"
          });
        }
      };

      for (let i = 0; i < numRepeats; i++) {
        for (let j = 0; j < numRepeats; j++) {
          for (let k = 0; k < numRepeats; k++) {
            drawUnitCell(i * a, j * a, k * a, a);
          }
        }
      }

      for (let i = 0; i <= numRepeats; i++) {
        for (let j = 0; j <= numRepeats; j++) {
          for (let k = 0; k <= numRepeats; k++) {
            viewer.addArrow({
              start: { x: i * a, y: j * a, z: k * a },
              end: { x: i * a, y: j * a, z: k * a + 1.5 },
              color: "red", radius: 0.3
            });
          }
        }
      }

      for (let i = 0; i < numRepeats; i++) {
        for (let j = 0; j < numRepeats; j++) {
          for (let k = 0; k < numRepeats; k++) {
            viewer.addArrow({
              start: { x: (i + 0.5) * a, y: (j + 0.5) * a, z: (k + 0.5) * a },
              end: { x: (i + 0.5) * a, y: (j + 0.5) * a, z: (k + 0.5) * a - 3.5 },
              color: "red", radius: 0.3
            });
          }
        }
      }

      viewer.zoomTo();
      viewer.spin(true);
      viewer.resize();
      viewer.render();
    }
    </script>


Initialize
----------

.. code-block:: python

   import magnon


Structure Inputs
----------------

Magnon makes good use of ASE, and this includes its I/O capabilities. The unit cell structure
can be read directly into an Atoms object using e.g. :code:`ase.io.read()` with any of the standard structure file formats allowed by that function.
In this example, we read from a POSCAR file :code:`FeGd.poscar`:

.. code-block::

    FeGd
    1
            1.0000000000         0.0000000000         0.0000000000
            0.0000000000         1.0000000000         0.0000000000
            0.0000000000         0.0000000000         1.0000000000
       Fe	Gd
        1	 1
    Direct
         0.000000000         0.000000000         0.000000000
         0.500000000         0.500000000         0.500000000

.. code-block::

   atoms = ase.io.read("FeGd.poscar")

Since not all of these file formats will allow the inclusion of magnetic moment information, and to allow flexibility in
the set-up of these quantities, the magnetic moments may be read in from a separate file using :code:`magnon.input.read_site_spin_data`,
which takes a file containing the magnitudes and directions of the magnetic moments, :code:`FeGd_moments`,

.. code-block:: text

   # magmom sx sy sz
     3.0     0  0  1
     7.0     0  0 -1

.. code-block::

   magmoms, directions = magnon.input.read_site_spin_data("FeGd_moments")
   atoms.set_initial_magnetic_moments(magmoms[:, np.newaxis] * directions)

The spin quantum number or magnetic moment may be specified by changing the *magmom_scaling* factor. To specify magnetic moments in
the input file, the scaling should be :math:`1`. When specifying spin quantum numbers, it should be :math:`2`.

Exchange coupling inputs
------------------------

The exchange coupling may be imported from a file :code:`FeGd_exchange` with the format,

.. code-block:: text

   # i j  vector        J
     1 1  1.0 0.0 0.0  18
     2 2  1.0 0.0 0.0   5
     1 2  0.5 0.5 0.5  -2

This format matches that returned by many codes used to compute exchange couplings. This is read in the form required to initialise an :code:`InteractionList` via
:code:`magnon.interactions.read_interactions`:

.. code-block::

   interactions = magnon.interactions.read_interactions("FeGd_exchange")
   interactions = magnon.interactions.InteractionList(interactions, atoms=atoms)

.. note::

   Magnon also supports tensor coupling for more complicated spin interactions. These are specified by nine entries after
   the vector, rather than the single entry for a scalar (isotropic) interaction.

Reading in a single line
------------------------

All of the above can be achieved more compactly using

.. code-block::

   atoms, interactions = magnon.input.create_interacting_system('FeGd.poscar', 'FeGd_moments', 'FeGd_exchange',)

Getting the bandstructure
-------------------------

The rest of the process proceeds along the lines we have already seen in previous tutorials; here we manually set the high-symmetry points as
we did in the very first tutorial:

.. code-block::

    interactions = interactions.symmetrize(atoms)

    special = {
        'G': [0,0,0],
        'H': [0.5,-0.5,0.5],
        'N': [0,0,0.5],
        'P': [0.25,0.25,0.25],
    }

    path = atoms.get_cell().bandpath(path='GHNGPH', npoints=240, special_points=special)

    spectrum = magnon.MagnonSpectrum(atoms, interactions, ham_prefactor=1)

    bstruct = spectrum.get_band_structure(path)
    bstruct.plot(emin=0, emax=200, filename='FeGd_bands.png')

yielding

.. image:: FeGd_bandstructure.png

Extracting the bandstructure data
---------------------------------

We can extract the reciprocal space path points and corresponding band energies from the BandStructure object as:

.. code-block::

   k_points = bstruct.path.cartesian_kpts()
   bands = bstruct.energies

and we can write the bandstructure to a JSON file using:

.. code-block::

   bstruct.write("FeGd_bands.json")

which allows you to use the raw data to generate your own plots or analyses.

Full script
-----------

.. code-block::

    import magnon

    atoms, interactions = magnon.input.create_interacting_system('FeGd.poscar', 'FeGd_moments', 'FeGd_exchange',)
    interactions = interactions.symmetrize(atoms)

    special = {
        'G': [0,0,0],
        'H': [0.5,-0.5,0.5],
        'N': [0,0,0.5],
        'P': [0.25,0.25,0.25],
    }

    path = atoms.get_cell().bandpath(path='GHNGPH', npoints=240, special_points=special)

    spectrum = magnon.MagnonSpectrum(atoms, interactions, ham_prefactor=1)

    bstruct = spectrum.get_band_structure(path)
    bstruct.plot(emin=0, emax=200, filename='FeGd_bands.png')

    k_points = bstruct.path.cartesian_kpts()
    bands = bstruct.energies

    print(k_points[0], bands[0][0])

    bstruct.write("FeGd_bands.json")

