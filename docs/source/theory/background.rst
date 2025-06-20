Background
===========

In this section, we cover the background theory needed to appreciate the derivation of the Linear Spin Wave Hamiltonian in the
next section. This is at the level of an undergraduate course in introductory quantum mechanics, and can serve as background reading
or a refresh!

We will cover:

* orbital and spin angular momentum of electrons
* the connection between angular momentum and magnetic moments
* quantum mechanical spin operators and their properties
* symmetries of multiple spin states
* the origins of the Heisenberg Hamiltonian

Introduction to magnetic moments
--------------------------------

Magnetic atoms are those with nonzero **magnetic moment**; the magnetic moment of an atom determines the strength of the
magnetic field it produces, and how strongly it will interact with other magnetic fields.

The magnetic moment is a vector property, i.e. it points in a certain direction. Often you will see magnetic moments
represented as arrows drawn on an atom, pointing in the same direction as the magnetic moment.

We start by examining the physical origins of the magnetic moments of atoms.

Electron orbits
++++++++++++++++

There are two sources of magnetic moment in atoms. The first is the **orbital angular momentum** of the electrons. This is momentum associated
with moving in an orbit, rather than a straight line. In this case, the size of the magnetic moment is determined by two things :cite:`griffiths_EM`:

1. The rate at which the electron orbits the atom, determined by the electron current.
2. The size of the orbit, described by the radius or the area traced out.


.. figure:: classical_moment.png
   :scale: 30 %

   *A magnetic moment represented as a vector whose magnitude depends on the area enclosed by an orbiting charge.*


For an electron orbiting an atom once every :math:`T` seconds at a distance :math:`r` with velocity :math:`v`, the size of the current,

.. math::
   :label: electron_orbits_current

   I = \frac{e}{T} = \frac{ev}{2\pi r}

with :math:`e` the electron charge. The area enclosed is :math:`\pi r^2`, so the size of the magnetic moment :math:`\mu` due to the orbital motion is

.. math::
   :label: electron_orbits_moment

   \mu = \frac{evr}{2}.

The orbital angular momentum is defined as :math:`L = mvr`, so

.. math::
   :label: electron_orbits_moment_L

   \mu = \frac{e}{2m}L.

In quantum mechanics, angular momentum is quantized in units of :math:`\hbar`, so we define a useful unit, the **Bohr magneton**,

.. math::
   :label: bohr_magneton

   \mu_B = -\frac{e\hbar}{2m},

and then the magnetic moment is some multiple of this unit,

.. math::
   :label: electron_orbit_mu_bohr_magneton

   \mu = \tilde{L}\mu_B,

where :math:`\tilde{L}` is the dimensionless angular momentum quantum number. In most references, the tilde is omitted.

We see that when an electron orbits an atom, its orbital motion generates a magnetic field which we quantify in terms of
a magnetic moment. However, there are other sources of angular momentum within the atom which are more fundamental.

Spin
++++

Experiments have shown that electrons have additional angular momentum which is not associated with their orbital
motion. In fact, this is intrinsic angular momentum, not associated with any type of motion, and it is called **spin angular momentum** :cite:`griffiths_QM`.

.. note::

   The name 'spin' was established by early proposals which suggested that the newly-discovered angular momentum came from the
   electron spinning about some axis, like a spinning top. This is now known not to be the case, as it would require the
   surface of the electron to move faster than the speed of light!

The spin angular momentum leads to a magnetic moment also, with a similar (but modified) form,

.. math::
   :label: spin_moment

   \mu = gS\mu_B,

where :math:`S` is the spin angular momentum. :math:`g` is the **spin g-factor**, which is approximately :math:`2`. The
appearance of this additional factor is not too surprising; our original choice of unit, the Bohr magneton, was derived by
considering orbiting electrons - but the spin angular momentum has nothing to do with orbits, so there is no reason
why it should be similarly expressed as a straightforward multiple only of :math:`\mu_B`.

.. note::

   It is common to refer to either a 'magnetic moment' or a 'spin' when discussing the magnetic properties of atoms. Both are vectors pointing in the same direction - the only difference is their magnitude, and physical units; spin is an angular momentum, not a magnetic moment.

We are now familiar with magnetic moments, and in Magnos we focus on the moment due to spin angular momentum only. But
to model the spin angular momentum accurately, we need to switch to the quantum mechanical formalism for magnetic moments.

Spin in quantum mechanics
--------------------------

We now associate each component of the spin angular momentum with a quantum mechanical operator, denoted by
:math:`\hat{S}^{x}`, :math:`\hat{S}^{y}` and :math:`\hat{S}^{z}`.

Commutator
++++++++++++

The most fundamental property of spin angular momentum operators in quantum mechanics is their **commutators**,

.. math::
   :label: spin_commutators

   \begin{aligned}
       \left[\hat{S}^x, \hat{S}^y\right] &= i\hbar \hat{S}^z \\
       \left[\hat{S}^y, \hat{S}^z\right] &= i\hbar \hat{S}^x \\
       \left[\hat{S}^z, \hat{S}^x\right] &= i\hbar \hat{S}^y
   \end{aligned}.


These commutators can be derived in two ways:

1. Considering that the spin angular momentum is a type of angular momentum operator, it must have the same commutation relations as orbital angular momentum by analogy. These can be derived using :math:`\hat{\mathbf{L}} = \hat{\mathbf{r}} \times \hat{\mathbf{p}}`
2. The most fundamental definition of angular momentum is as the **generator** of rotations in 3D space (meaning any rotation can be built up from infinitesimal rotations, which are described by the angular momentum operators). Then the commutator follows directly from the non-commutativity of 3D rotations around different axes in 3D space.

These commutation relations mean that two or more components of the spin angular momentum are not simultaneously well-defined.

Total angular momentum operator
+++++++++++++++++++++++++++++++

We also define the operator

.. math::
   :label: spin_total_operator

   \hat{\mathbf{S}}^2 = (\hat{S}^x)^2 + (\hat{S}^y)^2 + (\hat{S}^z)^2,

and it is easily shown that this commutes with each of the operators for the x,y and z components.

Ladder operators
++++++++++++++++

Define the ladder operators

.. math::
   :label: spin_ladder_operators

   \hat{S}^{\pm} = \hat{S}^x \pm i\hat{S}^y.

We will use this to understand the eigenstates of :math:`\hat{\mathbf{S}}^2` and :math:`\hat{S}^z` and their corresponding
quantum numbers.

Quantum number #1: :math:`m`
++++++++++++++++++++++++++++

The quantum number :math:`m` corresponds to the measurement of a component of the spin angular momentum along some fixed axis,
henceforth taken as the z-axis. We consider a state :math:`\ket{m}` which is an eigenstate of :math:`\hat{S}^z` with eigenvalue :math:`m\hbar`. Using the commutators,

.. math::
   :label: spin_m_number

   \hat{S}^{z} \hat{S}^{\pm} \ket{m} = (m\pm 1) \hat{S}^{\pm} \ket{m},

so that the states :math:`\hat{S}^{\pm} \ket{m}` may be identified with new states :math:`\ket{m\pm1}` and it is clear that
the operators :math:`\hat{S}^{\pm}` move us up or down the 'ladder' of eigenstates of :math:`\hat{S}^z`.

Quantum number #2: :math:`s`
+++++++++++++++++++++++++++++

The spin angular momentum is finite, so there must be an upper and lower bound on the eigenvalues the z-component of spin can adopt.
We'll therefore consider the state :math:`\ket{s}`, where :math:`s\hbar` is the maximum z-component of spin angular
momentum.

Because this is the maximum value, applying :math:`S^+` should destroy this state and give zero. Using the relations above,

.. math::
   :label: spin_s_plus_s_minus

   \hat{S}^-\hat{S}^+ = \hat{\mathbf{S}}^2 - (\hat{S}^z)^2 + \hbar\hat{S}^z.

Applying this to :math:`\ket{s}`, the left-hand side is zero because this is the maximum value allowed and further application
of :math:`S^+` destroys the state. Then

.. math::
   :label: spin_s_quantum_number

   \hat{\mathbf{S}}^2\ket{s} = \hbar s(s+1)\ket{s}.

Eigenspectrum
++++++++++++++

.. figure:: quantised_spin.png
   :figwidth: 400

   *A visualisation of spin states with fixed spin magnitude and determinate z-component, but indeterminate x- and y-components. The z-component is quantised in units of :math:`\hbar`.*

The spin angular momentum states are defined by the quantum numbers :math:`s` and :math:`m`:

* :math:`s` is fixed by the species of particle. For example, electrons have :math:`s = \frac{1}{2}`.
* :math:`m` determines the z-component of the spin angular momentum, and may take values from :math:`-s` to :math:`+s`

.. note::

   Notice that the z-component can never be equal to the length of the spin angular momentum vector (:math:`S^z_{max} = s\hbar < s(s+1)\hbar`) under the quantisation.
   This is necessary to satisfy the uncertainty relations above. If it could be equal, all components would be deterministic; :math:`S^z` would be equal to the spin
   of the particle, the other components would be zero.

Multiple spins
++++++++++++++

Consider a state of two spins, which have only :math:`m=\pm\frac{1}{2}`; there are then four unique states:

.. math::
   :label: two_spin_states

   \ket{\uparrow\uparrow}, \ket{\uparrow\downarrow}, \ket{\downarrow\uparrow}, \ket{\downarrow\downarrow},

and we wish to categorise these using the quantum numbers :math:`s` and :math:`m`. Components of angular momentum along
a certain axis is additive, so e.g.

.. math::
   :label: additive_spin_operators

   \hat{S}^z_{tot} = \hat{S}^z_{1} + \hat{S}^z_{2}.

So we can categorise the four states under :math:`m` easily enough.

In quantum mechanics, particles such as electrons are indistinguishable. This means that if we have two electrons, it is not
possible to say which one is which. This symmetry must be reflected in the quantum mechanical wavefunction when we swap the order of the electrons.
The states :math:`\ket{\uparrow\uparrow}` and :math:`\ket{\downarrow\downarrow}` satisfy this already, but the other two must be replaced by
symmetric linear combinations,

.. math::
   :label: two_spin_states_symmetric

   \frac{1}{\sqrt{2}}\left(\ket{\uparrow\downarrow} + \ket{\downarrow\uparrow}\right)

   \frac{1}{\sqrt{2}}\left(\ket{\uparrow\downarrow} - \ket{\downarrow\uparrow}\right).

When assigning :math:`s` quantum numbers, it is clear that :math:`\ket{\uparrow\uparrow}` and :math:`\ket{\downarrow\downarrow}` must correspond to
:math:`s = 1` since :math:`m = 1`, and :math:`s > m`. But understanding which of the above two corresponds to :math:`s=0,m = 0` and which one corresponds
to :math:`s=1,m = 0` is more difficult. Rather than applying :math:`\hat{\mathbf{S}}^2` directly, it is easier to consider
that the :math:`\hat{S}^{\pm}` do not change the :math:`s` quantum number, so by applying one of these to either :math:`\ket{\uparrow\uparrow}` or
:math:`\ket{\downarrow\downarrow}` we can find the other states with the same :math:`s` quantum number. In the two-spin example,
it turns out that the :math:`s=1` states are

.. math::
   :label: two_spin_states_s_1

   \ket{\uparrow\uparrow}

   \frac{1}{\sqrt{2}}\left(\ket{\uparrow\downarrow} + \ket{\downarrow\uparrow}\right)

   \ket{\downarrow\downarrow},

leaving

.. math::
   :label: two_spin_states_s_0

   \frac{1}{\sqrt{2}}\left(\ket{\uparrow\downarrow} - \ket{\downarrow\uparrow}\right)

as the only state with :math:`s=0`.

The significant result is that, under **exchange** of spins (i.e. swapping the positions of the first and second arrow),
the first set of states are symmetric (they don't change) and the second set (of one state) is antisymmetric (it is multiplied
by :math:`-1`). This means that there is a clear link between the **symmetry of the spin state** and the **total spin**. We will exploit this
property to great effect in setting up the Heisenberg Hamiltonian.

.. note::

   To see the above result for the more general case, consider that for multiple spins the operators :math:`\hat{S}^x`, :math:`\hat{S}^y`, :math:`\hat{S}^z` are additive
   i.e. :math:`\hat{S}^x = \hat{S}^x_1 + \hat{S}^x_2 + ...`, meaning the :math:`\hat{S}^{\pm}` are additive also. Thus application of
   :math:`\hat{S}^{\pm}` to the top or bottom state to get those inbetween (with the same :math:`s` quantum number) will preserve
   the symmetry or antisymmetry.

Exchange interaction
--------------------

Origins in symmetry
++++++++++++++++++++

The interaction we consider in Magnos is the **exchange interaction**. This is a purely quantum phenomenon, without a classical equivalent.

The indistinguishability of electrons requires that the spatial part of the wavefunction (that which depends on particle
positions) must show symmetry when they are exchanged. Consider two electrons in two quantum states labelled :math:`\psi_A` and :math:`\psi_B`. We cannot write the combined
wavefunction as

.. math::
   :label: exchange_symmetry_violation

   \psi_A(\mathbf{r}_1)\psi_B(\mathbf{r}_2),

because this suggests that we know which particle is described by each state :math:`A,B`. Instead, we use a combination of
terms in which we *exchange* the coordinates,

.. math::
   :label: exchange_symmetry_respected

   \psi_A(\mathbf{r}_1)\psi_B(\mathbf{r}_2) \pm \psi_A(\mathbf{r}_2)\psi_B(\mathbf{r}_1),

(ignoring for now normalisation). You can see that there are two options:

* '+' gives a *symmetric* joint wavefunction under :math:`\mathbf{r}_1\leftrightarrow\mathbf{r}_2`, because the new wavefunction is the same
* '-' gives an *antisymmetric* wavefunction under :math:`\mathbf{r}_1\leftrightarrow\mathbf{r}_2` because the new wavefunction is a negative copy of the original.

A system of electrons, which are fermions, must have a wavefunction which is antisymmetric under exchange overall. The overall wavefunction
is a product of one of the joint spatial wavefunctions above and the joint spin wavefunction, describing the orientations
of the spins. The spatial wavefunction chooses the lowest energy (either the symmetric
or antisymmetric) option and then the spin wavefunction adopts the opposite, to give overall antisymmetry.

Now consider taking, for a pair of spins,

.. math::
   :label: two_spins_total_spin

   (\mathbf{S}_1 + \mathbf{S}_2)^2 = (\mathbf{S}_1)^2 + (\mathbf{S}_2)^2 + 2 \mathbf{S}_1 \cdot \mathbf{S}_2,

which means

.. math::
   :label: two_spins_dot_product

   \mathbf{S}_1 \cdot \mathbf{S}_2 = \frac{1}{2} \left[ (\mathbf{S}_{tot})^2 - (\mathbf{S}_1)^2 - (\mathbf{S}_2)^2 \right],

from which we see that the dot product of two spin vectors depends on the :math:`s` quantum number of the joint spin state. (this is the eigenvalue of :math:`(\mathbf{S}_{tot})^2`)
However, we have already showed that:

* this quantum number determines the symmetry of the spin part of the wavefunction (under the exchange of spins),
* the symmetry of the spin part is determined by the symmetry of the spatial part, to give overall antisymmetry,
* the symmetry of the spatial part depends on which electron state is of lower energy.

**So there is a direct correlation between the dot product of two spin vectors and the energy of multi-electron
wavefunctions.** This underpins the rest of what we will do - the exchange interaction between spins is actually a way of
representing the energy differences between different spatial electron wavefunctions which we 'select' by playing with the directions
of the spin vectors :cite:`blundell2021`.

Heisenberg Hamiltonian
++++++++++++++++++++++

Using the above, we consider pairs of spins and quantify the effect on the energy of changing their relative orientation by the **exchange
interaction parameter** :math:`J`,

.. math::
   :label: basic_heisenberg_hamiltonian

   H = \sum_{\mathbf{R},\mathbf{R'}}\sum_{bb'} J_{\mathbf{R}b\mathbf{R'}b'} \mathbf{S}_{\mathbf{R}b} \cdot \mathbf{S}_{\mathbf{R'}b'},

known as the **Heisenberg Hamiltonian**. The sum runs over all unique pairs of atoms; :math:`b,b'` are indices of the atoms
in the unit cell, and :math:`\mathbf{R},\mathbf{R'}` are indices over cells (labelled by lattice translation vectors). :math:`J` is the exchange coupling interaction strength,
and taking the dot product between the :math:`\mathbf{S}` quantifies the relative orientation of pairs of spin
vectors.

