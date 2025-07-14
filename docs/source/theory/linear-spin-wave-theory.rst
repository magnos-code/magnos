Linear Spin Wave Theory
========================

Heisenberg Hamiltonian
----------------------

General Heisenberg Hamiltonian
+++++++++++++++++++++++++++++++

The starting point to describe magnon excitations is the Heisenberg Hamiltonian written in general form (applicable to any material), 

.. math::
   :label: general_heisenberg_hamiltonian

   H = -\frac{1}{2}\sum_{\mathbf{R},\mathbf{R'}} \sum_{b,b'} \sum_{\alpha, \alpha'} \hat{S}^{\alpha}_{\mathbf{R}b} J^{\alpha\alpha'}_{\mathbf{R}n\mathbf{R'}b'} \hat{S}^{\alpha'}_{\mathbf{R'}b'},

where :math:`\mathbf{R}, \mathbf{R'}` are direct lattice vectors, :math:`b,b'` are indices that label magnetic sites in the unit cell, and
:math:`\alpha,\alpha'` are the Cartesian indices :math:`{x,y,x}`. :math:`J` is the strength of the exchange coupling interaction
in units of energy, and :math:`\hat{S}^{\alpha}` is the quantum mechanical operator associated with the :math:`\alpha`-component
of spin angular momentum, in units of :math:`\hbar`.

.. figure:: R_tau.png
   :figwidth: 400

   *Explanation of the notation used. The position vector of a site is split into a lattice vector part and a unit cell part.*

The inclusion of the :math:`\alpha,\alpha'` indices allows the interaction between spins to be anisotropic and include off-diagonal
components. The simplest case is an identity matrix in these indices, multiplied by the magnitude of the interaction; this corresponds
to the basic form given in the previous section. The more general form allows us to include:

**Direct exchange**
To maximize the separation of electrons (which repel under the Coulomb interaction) the antisymmetric spatial wavefunction is preferred.
This therefore favours a symmetric spin wavefunction, with spins aligned along the same direction (ferromagnetism).

This interaction populate the elements :math:`J_{\mathbf{R}b\mathbf{R}'b'}^{\alpha\alpha}`.

**Indirect exchange**
The electrons of some intermediate part of the system affect
the energies of the two spin states. For example, the superexchange mechanism, which is responsible for antiferromagnetism
in some materials, involves electron interaction via an intermediate nonmagnetic atom.

The interaction terms are the same as for direct exchange, where :math:`\alpha=\alpha'`

**Kitaev interaction**
An anisotropic, bond-directional interaction on 2D honeycomb lattices.

Depending on the direction of the bond, the interaction either takes the form :math:`J_{\mathbf{R}b\mathbf{R'}b'}^{xx}S_{\mathbf{R}b}^x S_{\mathbf{R'}b'}^x`,
:math:`J_{\mathbf{R}b\mathbf{R'}b'}^{yy}S_{\mathbf{R}b}^y S_{\mathbf{R'}b}^y` or :math:`J_{\mathbf{R}b\mathbf{R'}b'}^{zz}S_{\mathbf{R}b}^z S_{\mathbf{R'}b'}^z`


**Dzyaloshinskii-Moriya interaction**
An interaction associated with the vector product :math:`\mathbf{S}_1 \times \mathbf{S}_2` of two spin vectors, resulting
from the spin-orbit interaction.

This interaction populate the elements :math:`J_{\mathbf{R}b\mathbf{R}'b'}^{\alpha\alpha'}` which are off-diagonal i.e. :math:`\alpha\ne\alpha'`

Holstein-Primakoff Transformation
----------------------------------

To describe quantized magnetic excitations, we want to write this Hamiltonian in terms of **creation and annihilation
operators** in the formalism of **second quantization**.

Creation and annihilation operators
++++++++++++++++++++++++++++++++++++

.. figure:: spin_flip.png
   :figwidth: 600

   Figure 1. The effect of the creation operator :math:`\hat{b}^{\dagger}_{\mathbf{R}b}` on a ferromagnetic ground state.

A magnetic excitation causes deviation from a perfect ferromagnetic arrangement where all spins point in the same direction.
Figure 1 shows the effect of a spin-flip operator on such a ground state.

* The operator :math:`\hat{b}^{\dagger}_{\mathbf{R}b}` creates a magnetic excitation at the site located at :math:`\mathbf{R}+\boldsymbol{\tau}_b`.
* The operator :math:`\hat{b}_{\mathbf{R}b}` destroys a magnetic excitation at the site located at :math:`\mathbf{R}+\boldsymbol{\tau}_b`.

The creation and annihilation operators are akin to ladder operators, which we have already seen in the discussion of the
spin angular momentum. Here, the creation and annihilation operators act on **number states**, and raise or lower the number of particles by one by
moving between adjacent number states. A number state has a fixed number
of particles; in this case a 'particle' is a magnetic excitation (more accurately a '**quasiparticle**')

Magnetic excitations change the spin by :math:`\hbar`, so they are quasiparticles of integer spin, which makes them bosons (fermions have half-integer spin).
The important property of these operators in order that they describe **bosons** is that they must have the commutation relation

.. math::
   :label: boson_commutation_relations

   \left[\hat{b}_{\mathbf{R}b}, \hat{b}^{\dagger}_{\mathbf{R'}b'} \right] = \delta_{\mathbf{R}\mathbf{R'}}\delta_{bb'}.

From this, we may show directly by action on a **number state** that the operator :math:`\hat{n}_{\mathbf{R}b} = \hat{b}^{\dagger}_{\mathbf{R}b}\hat{b}_{\mathbf{R}b}`
can be used to find the number of magnetic excitations present (the **number operator**).

Magnetic excitations
++++++++++++++++++++++

The ferromagnetic ground state corresponds to maximum z-component of spin, so a magnetic excitation will lower this component,
reducing the :math:`m` quantum number.
We may then be tempted to associate the spin raising and lowering operators with the creation and annihilation operators above as

.. math::
   :label: naive_operators

   \hat{S}^{+}_{\mathbf{R}b} \sim \hat{b}_{\mathbf{R}b}

   \hat{S}^{-}_{\mathbf{R}b} \sim \hat{b}^{\dagger}_{\mathbf{R}b},

but doing so results in the wrong commutator for the spin component operators :math:`\hat{S}^x`, :math:`\hat{S}^y` and :math:`\hat{S}^z`.

Instead, we must use the transformation due to Holstein and Primakoff :cite:`holstein1940`,

.. math::
   :label: holstein_primakoff

   \hat{S}^{+}_{\mathbf{R}b} = \sqrt{2 S_{\mathbf{R}b}} \left( 1 - \frac{\hat{b}^{\dagger}_{\mathbf{R}b}\hat{b}_{\mathbf{R}b}}{2 S_{\mathbf{R}b}} \right)^{\frac{1}{2}} \hat{b}_{\mathbf{R}b}

   \hat{S}^{-}_{\mathbf{R}b} = \sqrt{2 S_{\mathbf{R}b}} \hat{b}^{\dagger}_{\mathbf{R}b} \left( 1 - \frac{\hat{b}^{\dagger}_{\mathbf{R}b}\hat{b}_{\mathbf{R}b}}{2 S_{\mathbf{R}b}} \right)^{\frac{1}{2}}

which, it can be shown, results in the correct commutation relations for the components of spin. In addition, we see that the form
of the bracket enforces that there is a minimum :math:`m` state as it will be zero when we try to excite a number of magnetic quasiparticles
which exceeds :math:`2 S_{\mathbf{R}b}`. The z-component of spin will be

.. math::
   :label: holstein_primakoff_z

   \hat{S}^z_{\mathbf{R}b} = S_{\mathbf{R}b} - \hat{b}^{\dagger}_{\mathbf{R}b}\hat{b}_{\mathbf{R}b} = S_{\mathbf{R}b} - \hat{n}_{\mathbf{R}b}.

If the perturbation to the ground state is small, so that the expected number of magnetic excitations in a state satisfies

.. math::
   :label: holstein_primakoff_approx_condition

   \frac{n_{\mathbf{R}b}}{2 S_{\mathbf{R}b}} \ll 1,

then we may use the binomial approximation on the bracket to obtain

.. math::
   :label: holstein_primakoff_approx

   \hat{S}^{+}_{\mathbf{R}b} = \sqrt{2 S_{\mathbf{R}b}}  \hat{b}_{\mathbf{R}b}

   \hat{S}^{-}_{\mathbf{R}b} = \sqrt{2 S_{\mathbf{R}b}} \hat{b}^{\dagger}_{\mathbf{R}b},

which is nothing more than our naive guess with an extra prefactor! In this approximation we still break the commutation
relations for the spin components, but the perturbation to the ground state is small enough for this to have a negligible effect.
This is analogous to considering the commutation of rotations in 3D space; the smaller the size of the rotations, the less important
the order in which they are applied will be.
Therefore, although we have obtained the same result as our original guess, we better understand the conditions on its use.

Given that spin is quantized, you may (very justifiably) feel uncomfortable with this last approximation. Since there is a minimum change of spin of :math:`\hbar`
associated with a magnetic excitation, there appears to be a restriction preventing the perturbation from being arbitrarily small.
However, we will later find that the eigenstate excitations are not localized spin excitations, but Fourier transformed excitations.
This corresponds to an excitation smeared out over many sites, with a small average deviation on any individual site.


Rotated frame approach
-----------------------

We have obtained a transformation from the spin component operators to creation and annihilation operators
which excite a single spin away from the ferromagnetic ground state. However, we wish to consider
arbitrary spin ground states in which the spins are not necessarily aligned in the same direction. To do this, we use the
rotated frame approach due to Toth and Lake :cite:`toth2015`.

.. figure:: rotated_frame.png

   *Figure 2. A schematic showing the two types of rotation used to describe different ground sate spin configurations.*

We retain the operators we saw in the last section, and map the more general ground state onto a ferromagnetic ground state
to allow us to use these operators. This requires expressing the actual ground state as a ferromagnetic ground state in which
the spins are rotated by two types of rotation:

* A **site rotation**, :math:`Q^{site}_{b}` which determines how each spin in the unit cell at the origin is rotated from the ferromagnetic state.
* A **cell rotation**, :math:`Q^{cell}_{\mathbf{R}}` which determines how much all of the spins are rotated between adjacent lattice cells.

Figure 2 shows these two rotations in action. We then rewrite the Heisenberg Hamiltonian using these rotations,

.. math::
   :label: rotated_spin

   \mathbf{S}_{\mathbf{R}b} = Q^{cell}_{\mathbf{R}}Q^{site}_{b}S_{\mathbf{R}b}^{FM},

where :math:`S_{\mathbf{R}b}^{FM}` are the spin vectors of the ferromagnetic ground state, all pointing along :math:`\hat{\mathbf{z}}`.

Site rotation
++++++++++++++

The site rotation can be expressed in column-vector form as

.. math::
   :label: site_rotation

   Q^{site}_{b} = [\mathbf{q}_{b,1}, \mathbf{q}_{b,2}, \mathbf{q}_{b,3}].

Since the ferromagnetic spin vector is just :math:`[0, 0, 1]` this means that :math:`\mathbf{q}_b^3` should be a unit vector
in the direction of the spin, and with the other two vectors should form an orthonormal basis.

To write the Hamiltonian in terms of :math:`\hat{b}_{\mathbf{R}b}` and :math:`\hat{b}_{\mathbf{R}b}^\dagger` we write in terms of
:math:`S_{\mathbf{R}b}^{\pm}` and :math:`S_{\mathbf{R}b}^{z}` and use the transformation discussed in the previous section. We find
that the factors of complex :math:`i` introduced by the conversion to :math:`S_{\mathbf{R}b}^{\pm}` makes it more natural to instead
use the vectors

.. math::
   :label: u_v_vectors

   \mathbf{u}_{b} = \mathbf{q}_{b,1} + i \mathbf{q}_{b,2}

   \mathbf{v}_{b} = \mathbf{q}_{b,3}.

These vectors are generated internally by Magnos from the spin direction vectors.

Cell rotation
+++++++++++++++

The cell rotation, being commensurate with the lattice, may be represented by a wavevector where the angle of rotation
in a given cell is

.. math::
   :label: cell_rotation_angle

   \theta = \mathbf{k}_{rot}\cdot\mathbf{R}.

A rotation axis should also be specified. Magnos takes :math:`\mathbf{k}_rot` and applies the magnetic ordering.

.. figure:: mag_ord.png

   *Figure 3. Two different ways to represent an antiferromagnet: left, using a magnetic ordering wavevector; right, using a supercell. The left approach is preferrable.*


If it is possible to represent a structure using a magnetic ordering wavevector, it is better to represent it this way than
manually building the ordering by inputting a supercell (the two options for a simple example are shown in Figure 3).
The supercell will not be the primitive cell, which will lead to spurious additional bands.

Eigensolution
--------------

In the previous section, we saw how to set up the Hamiltonian using creation and annihilation operators for different ground state
spin configurations. Now, we are ready to diagonalize. This will give us the eigenstates and eigenenergies of the Hamiltonian,
which are the longest-lived excitations under the model. This is because the Hamiltonian determines the time evolution of the system; an
eigenstate will be unchanged under time evolution.

Magnons
--------

As with many excitations, we find that we need to Fourier transform to get the eigenstates,

.. math::
   :label: operator_fourier_transform

   \hat{b}_{\mathbf{R}b} = \frac{1}{\sqrt{N}} \sum_{\mathbf{k}} \hat{b}_{\mathbf{k}b} e^{i\mathbf{k}\cdot (\mathbf{R}+\boldsymbol{\tau_b})}

   \hat{b}_{\mathbf{R}b}^\dagger = \frac{1}{\sqrt{N}} \sum_{\mathbf{k}} \hat{b}_{\mathbf{k}b}^\dagger e^{-i\mathbf{k}\cdot (\mathbf{R}+\boldsymbol{\tau_b})}.

The excitations created and destroyed by these operators are then **spin waves** or **magnons** rather than spin-flip excitations.
A magnon represents a spin flip smeared out over many spins, as shown in Figure 4c.

.. figure:: spin_excitations.png

   *Figure 4. a) A ferromagnetic ground state; b) a spin flip excitation; c) a magnon excitation*

Spin Wave Hamiltonian
----------------------

Rearranging and ensuring that we retain the Hermitian property of the Hamiltonian, its final block form is

.. raw:: html

   <div class="scroll-math">

.. math::
   :label: bdg_hamiltonian

   H =
   \begin{pmatrix}
   \hat{b}_{\mathbf{k}b}^\dagger & \hat{b}_{-\mathbf{k}b}
   \end{pmatrix}
   \begin{pmatrix}
   B_{1bb'}&B_{2bb'}\\
   B_{3bb'}&B_{4bb'}\\
   \end{pmatrix}
   \begin{pmatrix}
   \hat{b}_{\mathbf{k}b'} \\ \hat{b}_{-\mathbf{k}b'}^\dagger
   \end{pmatrix}

.. raw:: html

   </div>

with

.. raw:: html

   <div class="scroll-math">

.. math::
   :label: bdg_block_1

   B_{1bb'}(\mathbf{k}) = -\frac{1}{2}\sum_{\alpha,\alpha'}\sum_{\mathbf{R''}}\Bigg\{\frac{\sqrt{S_b S_{b'}}}{2} e^{i\mathbf{k}\cdot\mathbf{R''}}e^{i\mathbf{k}\cdot(\boldsymbol{\tau}_{b'}-\boldsymbol{\tau}_b)}u^\alpha_b\tilde{J}^{\alpha\alpha'}_{bb'\mathbf{R''}}\overline{u}^{\alpha'}_{b'} \\ + \sum_{b''}-S_{b''}\delta_{bb'}v^\alpha_{b''}\tilde{J}^{\alpha\alpha'}_{b''b'\mathbf{R''}}v^{\alpha'}_{b'}\Bigg\}

.. raw:: html

   </div>

   <div class="scroll-math">

.. math::
   :label: bdg_block_2

   B_{2bb'}(\mathbf{k}) = -\frac{1}{2}\sum_{\alpha,\alpha'}\sum_{\mathbf{R''}}\Bigg\{\frac{\sqrt{S_b S_{b'}}}{2} e^{i\mathbf{k}\cdot\mathbf{R''}}e^{i\mathbf{k}\cdot(\boldsymbol{\tau}_{b'}-\boldsymbol{\tau}_b)}u^\alpha_b\tilde{J}^{\alpha\alpha'}_{bb'\mathbf{R''}}u^{\alpha'}_{b'}\Bigg\}

.. raw:: html 

   </div>

   <div class="scroll-math">

.. math::
   :label: bdg_block_3

   B_{3bb'}(\mathbf{k}) = -\frac{1}{2}\sum_{\alpha,\alpha'}\sum_{\mathbf{R''}}\Bigg\{\frac{\sqrt{S_b S_{b'}}}{2} e^{i\mathbf{k}\cdot\mathbf{R''}}e^{i\mathbf{k}\cdot(\boldsymbol{\tau}_{b'}-\boldsymbol{\tau}_b)}\overline{u}^\alpha_b\tilde{J}^{\alpha\alpha'}_{bb'\mathbf{R''}}\overline{u}^{\alpha'}_{b'} \Bigg\}

.. raw:: html

   </div>

   <div class="scroll-math">

.. math::
   :label: bdg_block_4

   B_{4bb'}(\mathbf{k}) = -\frac{1}{2}\sum_{\alpha,\alpha'}\sum_{\mathbf{R''}}\Bigg\{\frac{\sqrt{S_b S_{b'}}}{2} e^{i\mathbf{k}\cdot\mathbf{R''}}e^{i\mathbf{k}\cdot(\boldsymbol{\tau}_{b'}-\boldsymbol{\tau}_b)}\overline{u}^\alpha_b\tilde{J}^{\alpha\alpha'}_{bb'\mathbf{R''}}u^{\alpha'}_{b'} \\ + \sum_{b''}-S_{b''}\delta_{bb'}v^\alpha_{b''}\tilde{J}^{\alpha\alpha'}_{b''b'\mathbf{R''}}v^{\alpha'}_{b'}\Bigg\}.

.. raw:: html
   
   </div>

This Hamiltonian is implemented in the MagnonSpectrum class in Magnos.

Bogoliubov transformation
-------------------------

Notation
+++++++++

Let's denote the above Hamiltonian

.. math::
   :label: bdg_compact

   H(\mathbf{k}) =  b^\dagger_{\mathbf{k}\rho}\mathcal{B}_{\rho\rho'}(\mathbf{k})b_{\mathbf{k}\rho'},

where

.. math::
   :label: bdg_compact_operator

   b_{\mathbf{k}\rho} =
   \begin{pmatrix}
   \hat{b}_{\mathbf{k}b'} \\ \hat{b}_{-\mathbf{k}b'}^\dagger
   \end{pmatrix}

and :math:`\mathcal{B}` is the 'Bogoliubov-de-Gennes'-type matrix,

.. math::
   :class: no-scroll
   :label: bdg_compact_matrix

   \mathcal{B} =
   \begin{pmatrix}
   B_{1bb'}&B_{2bb'}\\
   B_{3bb'}&B_{4bb'}\\
   \end{pmatrix}.

Original commutator
+++++++++++++++++++

We now check the commutator,

.. math::
   :label: original_commutator

   [b_{\mathbf{k}\rho},b_{\mathbf{k}\rho'}^\dagger]
   = \begin{pmatrix}
   [\hat{b}_{\mathbf{k}b},\hat{b}^\dagger_{\mathbf{k}b'}]&[\hat{b}_{\mathbf{k}b},\hat{b}_{-\mathbf{k}b'}] \\
   [\hat{b}^\dagger_{-\mathbf{k}b},\hat{b}^\dagger_{\mathbf{k}b'}]&[\hat{b}^\dagger_{-\mathbf{k}b},\hat{b}_{-\mathbf{k}b'}]\\
   \end{pmatrix}

   = \begin{pmatrix}
   \delta_{bb'}&0\\
   0&-\delta_{bb'}\\
   \end{pmatrix}

   = \tilde{\delta}_{\rho\rho'},

where :math:`\delta` is the identity matrix and :math:`\tilde{\delta}` is the **paraidentity** matrix, following Colpa.

Transformation conditions
+++++++++++++++++++++++++

We wish to transform the creation and annihilation operator basis to one in which the Hamiltonian is diagonal,

.. math::
   :label: diagonalised_hamiltonian

   H(\mathbf{k}) =  b^\dagger_{\mathbf{k}\sigma}\mathcal{T}_{\mathbf{k}\sigma\gamma}^\dagger\Big((\mathcal{T}_{\mathbf{k}}^\dagger)^{-1}_{\gamma\rho}\mathcal{B}_{\rho\rho'}(\mathbf{k})\mathcal{T}_{\mathbf{k}\rho'\gamma'}^{-1}\Big)\mathcal{T}_{\mathbf{k}\gamma'\sigma'}b_{\mathbf{k}\sigma'}

   = \beta_{\mathbf{k}\mu}^\dagger\mathcal{E}_{\mu\mu'}(\mathbf{k})\beta_{\mathbf{k}\mu'}.

However, we must ensure that the commutator above is preserved,

.. math::
   :label: commutator_condition

   [\beta_{\mathbf{k}\mu},\beta_{\mathbf{k}\mu'}^\dagger] = \tilde{\delta}_{\mu\mu'}.

It can be shown that this imposes the condition

.. math::
   :label: paraunitary_transformation_condition

   \sum_{\sigma,\sigma'}\mathcal{T}_{\mathbf{k}\mu\sigma}\tilde{\delta}_{\sigma\sigma'}\mathcal{T}^\dagger_{\mathbf{k}\sigma'\mu'} = \tilde{\delta}_{\mu\mu'}.

Diagonalization procedure
++++++++++++++++++++++++++

The diagonalization procedure which sustains this condition is due to Colpa :cite:`colpa1978`. First, we Cholesky decompose our Bogoliubov-de-Gennes-like Hamiltonian,

.. math::
   :label: cholesky_decomposition

   \mathcal{B}(\mathbf{k})_{\rho\rho'} = \sum_{\sigma} \mathcal{C}^\dagger_{\rho\sigma} \mathcal{C}_{\sigma\rho'}.

The Hermitian requirement is satisfied by properly symmetrising the couplings and selecting the correct ground state spin directions.
In order for this decomposition to be possible, the eigenvalues must also be positive. However, in some places they may be zero, so we
add small quantity to the diagonal,

.. math::
   :label: diagonal_small_quantity

   \mathcal{B}_{\rho\rho'} \to \mathcal{B}_{\rho\rho'}+\varepsilon\delta_{\rho\rho'}.

This can be removed from the obtained eigenvalues since the identity commutes with the Hamiltonian, so it has the same eigenstates.
We then construct the matrix

.. math::
   :label: U_matrix

   \mathcal{U}=\mathcal{C}\mathcal{T}^{-1}\mathcal{E}^{-\frac{1}{2}},

which, from the definition of the transformation, can be shown to be unitary,

.. math::
   :label: U_matrix_unitary_property

   \mathcal{U}^\dagger\mathcal{U} = \delta.

The eigenvalues can then be found by diagonalising :math:`\mathcal{C} \tilde{\delta}\mathcal{C}^\dagger`,

.. math::
   :label: finding_eigenvalues

   \mathcal{U}^\dagger[\mathcal{C} \tilde{\delta}\mathcal{C}^\dagger]\mathcal{U} = \tilde{\delta}\mathcal{E},

where :math:`\mathcal{E}` has two copies of the eigenvalues on the diagonal.

Duplicate eigenvalues
++++++++++++++++++++++

The reason for obtaining two copies of the eigenvalues is our use of the rotated frame approach. To describe systems beyond ferromagnetic
order using creation and annihilation operators for a ferromagnetic system, we had to introduce additional cross terms like :math:`\hat{b}\hat{b}`, :math:`\hat{b}\hat{b}^\dagger` and :math:`\hat{b}^\dagger\hat{b}^\dagger`, which leads to the
Bogoliubov-de-Gennes form of the matrix. In the ferromagnetic case, we end up with a block diagonal matrix where the two blocks on the
diagonal are equal and the problem reduces to the standard Hamiltonian with only a :math:`\hat{b}^\dagger\hat{b} = \hat{n}` term.

By considering the block diagonal form of the transformation matrix, it is possible to show that there are always two copies of identical
eigenvalues.

Limitations of Linear Spin Wave Theory
--------------------------------------

There are known limitations of using Linear Spin Wave Theory:

* The derivation assumes a small excitation number. This means the temperature should be low enough that there are not too many magnons.
* The Hamiltonian may lead to spurious symmetry enhancement - see Gohlke (2023) :cite:`gohlke2023`


