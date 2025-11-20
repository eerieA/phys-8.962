<!-- TOC -->

- [Lecture 5](#lecture-5)
    - [Dust element energy density](#dust-element-energy-density)
    - [Energy density of a moving frame P t and N t](#energy-density-of-a-moving-frame-p-t-and-n-t)
    - [Relating T tt with the full tensor](#relating-t-tt-with-the-full-tensor)

<!-- /TOC -->

# Lecture 5

## Dust element energy density

$\rho_0 = m n_0$

- m is particle's rest mass
- $n_0$ is rest density of particles in this element

$\rho_0$ is energy because rest energy of particles in a "traditional" view is $E=mc^2$, and in relativity we usually have $c=1$. This is why it **looks like** energy density = mass density.

Furthermore, in relativity:

$$
E = \gamma m, \quad \gamma = \frac{1}{\sqrt{1-v^2}}
$$

If the particles are at rest $ (v=0) $, so $\gamma = 1$ and:

$$
E = m.
$$

No kinetic energy. No momentum. Just rest energy.

## Energy density of a moving frame P t and N t

Related to this:

... 
Go into frame moving with v relative to this (rest) frame  
Energy density $\equiv \rho = (\gamma m) (\gamma n_0) = \gamma^2 \rho_0$.  
...  
We assembled $\rho$ by combining energy - timelike component of $\vec{p}$, with number density - timelike component of the number vector.
So these are 2 timelike components of 2 4-vectors here.  
...  
So we can write $\rho = P^t N^t$.

**What is $\vec{p}$ here?**

It is the **4-momentum *per particle***:

$
p^\mu = m u^\mu.
$

I guess the instructor mentally has it like $p^\mu = (P^t, P^i)$ too.

* The **timelike component** $P^t = P^0 = \gamma m$ is the **energy per particle** in the moving frame $in (c=1) units$.
* The spacelike components $P^i = \gamma m v^i$ are the particle’s momentum.

**What is $N^\mu$?**

This is a true 4-vector that represents:

* how many particles per unit volume, and
* how they are flowing.

$$
u^\mu_0 = (1,0,0,0)
\Rightarrow
N^\mu_0 = (n_0, 0,0,0).
$$

$$
N^\mu = ( \gamma n_0, \gamma n_0 v^i ).
$$

So $N^t$ is the boosted number density (because volumes Lorentz-contract).

**Why not use the true momentum density 4-vector?**

Because momentum *density* is **not** a 4-vector - densities transform with $\gamma$’s and are part of a rank-2 tensor (the stress-energy tensor).

**$\rho$ is not a scalar**

Because remember that in relativity, scalar means Lorentz invariant. And clearly $\rho$ depends on $\gamma$, i.e. frame dependent.

## Relating T tt with the full tensor

For motivation, we inspect the energy density in any frame:

$$
T^{tt} \propto u^t u^t.
$$

Notice the structure of it. There is a pattern we will discover.

**What about momentum density?**

In relativity, momentum density is the (not timelike) flux of energy, so it must be:

> But energy is timelike.  
> “timelike vs spacelike” only applies to 4-vectors, not rank-2 tensors.

$$
T^{ti}.
$$

Compute it from the two 4-vectors:

$$
T^{ti} = P^i N^t = (m u^i)(n_0 u^t)
= \rho_0 u^t u^i.
$$

Same structural pattern.

**What about momentum flux / stresses?**

A pressureless dust has no internal forces and no pressure. So the only flux of momentum (spatial) in direction $i$ across surface $j$ comes from the dust **carrying its own momentum along its own flow**.

Particle momentum:
$$
p^i = m u^i.
$$

Number flux (particles crossing surface per unit area):
$$
N^j = n_0 u^j.
$$

So:

$$
T^{ij} = P^i N^j
= (m u^i)(n_0 u^j)
= \rho_0 u^i u^j.
$$

Again the same pattern.

**Now the pattern is clear**

Every component of the stress–energy tensor for dust is built from:

* a factor of $m u^\mu$ (particle 4-momentum)
* a factor of $n_0 u^\nu$ (particle number current)

Multiplying gives:

$$
T^{\mu\nu} = (m n_0) u^\mu u^\nu = \rho_0 u^\mu u^\nu.
$$

This is exactly this outer product:

$$
T = N \otimes P.
$$

For each pair of indices $(\mu,\nu)$:

* $\nu$ tells you what kind of 4-momentum you're tracking (energy, momentum-x, etc.)
* $\mu$ tells you across which spacetime surface you're measuring the flow.

If particles:

* have 4-momentum $p^\nu$, and
* flow with 4-current $N^\mu$,

then the flux is:

$$
T^{\mu\nu} = N^\mu p^\nu.
$$

**In a nutshell**

Once we found $T^{tt} = \rho_0 (u^t)^2$, the rest follows automatically by covariance: A tensor is not allowed to have only one component defined in one frame.

If we know how $T^{tt}$ transforms, the entire tensor must be the outer product that produces that result.

Thus, the general covariant completion is:

$$
T^{\mu\nu} = \rho_0 u^\mu u^\nu.
$$

So both covariance and physical meaning require us to extend $T^{tt}$ into a full rank-2 object.