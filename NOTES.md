<!-- TOC -->

- [Lecture 1](#lecture-1)
    - [Deriving the Lorentz Transformations](#deriving-the-lorentz-transformations)
- [Lecture 2](#lecture-2)
    - [Kronecker delta and how to use for basis vecs](#kronecker-delta-and-how-to-use-for-basis-vecs)
    - [What does the eta (tensor) do in dot product](#what-does-the-eta-tensor-do-in-dot-product)
- [Lecture 3](#lecture-3)
    - [Contravariant covariant “upstairs” and “downstairs”](#contravariant-covariant-upstairs-and-downstairs)
    - [Simple contraction of tensors](#simple-contraction-of-tensors)
        - [References](#references)
        - [Even simpler 2 dim tensors example](#even-simpler-2-dim-tensors-example)
    - [Simple outer product example](#simple-outer-product-example)
    - [An index lowering example](#an-index-lowering-example)
        - [The Metric Tensor ($g{\mu\nu}$)](#the-metric-tensor-g\mu\nu)
        - [The Index Lowering Operation](#the-index-lowering-operation)
        - [Calculation of it](#calculation-of-it)
- [Lecture 4](#lecture-4)
    - [“Trajectory” of photon is null curve](#trajectory-of-photon-is-null-curve)
        - [Intuition and proper time](#intuition-and-proper-time)
        - [The physical reasoning](#the-physical-reasoning)
    - [Normalizing 4-velo?](#normalizing-4-velo)
    - [How mass connects to the norm](#how-mass-connects-to-the-norm)
    - [What is h omega](#what-is-h-omega)
    - [Levi–Civita and Permutation Parity](#levicivita-and-permutation-parity)
    - [Derive 4D Gauss Thm thru a V4 integral](#derive-4d-gauss-thm-thru-a-v4-integral)

<!-- /TOC -->

# Lecture 1

## Deriving the Lorentz Transformations

Tutorial on deriving it for $v$ along $x$ axis:  
https://www.youtube.com/watch?v=FvqutkaPmas

> I like the explaination on why the Lorentz Transformation was chosen to be linear. This instructor thought it is because of the principle of relativity:
>
> - The transformation should be "reversible".
> - Also the transformation from $K$ to $K'$ and its "reverse" should be of the same form (no "previledged" IRF needed to determine which is which), and a linear one satisfies that.

Tutorial on deriving it for $v$ along an arbitrary axis:  
https://www.youtube.com/watch?v=Afd34FuG65A


# Lecture 2

## Kronecker delta and how to use for basis vecs

Related to this part on the blackboard

$$
\begin{array}{rl}
     \vec{e}_{\alpha} & = \Lambda^{\bar{\beta}}_{\alpha} (\underset{\sim}{v}) \vec{e}_{\bar{\beta}} \\
        & = \Lambda^{\bar{\beta}}_{\alpha} (\underset{\sim}{v}) [ \Lambda^{\gamma}_{\bar{\beta}} (-\underset{\sim}{v}) \vec{e}_{\gamma}] \\
        & = [\Lambda^{\bar{\beta}}_{\alpha} (\underset{\sim}{v}) \Lambda^{\gamma}_{\bar{\beta}} (-\underset{\sim}{v})] \vec{e}_{\gamma}
\end{array}
$$

$$
\delta_{\alpha}^{\gamma} = \Lambda^{\bar{\beta}}_{\alpha} \Lambda^{\gamma}_{\bar{\beta}}
$$

We see a Kronecker delta $\delta_{\alpha}^{\gamma}$ appears at the end of that part.

So Kronecker delta is just defined as this

$$
\vec{e}_{\alpha} = \delta^{\gamma}_{\alpha} \vec{e}_{\gamma}
$$
, which roughly means

“The basis vector $\vec e_\alpha$ can be expressed as a linear combination of the same basis vectors ${\vec e_\gamma}$, with coefficients $\delta^\gamma{}_\alpha$."

It is $1$ if $\gamma=\alpha$, and $0$ otherwise. Like a "switch" or "binary picker".

Here is an example to get a more practical understanding of it. To write the sum of the trace of a 3 by 3 matrix using Kronecker delta:

$$
\text{Tr}A \doteq \sum_{i=1}^{3} \sum_{i=j}^{3} a_{ij}\delta_{ij}
= \sum_{i=1}^{3} a_{ii}
$$

. Following this it is easy to see that the identity matrix can be written as:

$$
I\doteq \begin{pmatrix}
\delta_{11} & \delta_{12} & \delta_{13}\\
\delta_{21} & \delta_{22} & \delta_{23}\\
\delta_{31} & \delta_{32} & \delta_{33}
\end{pmatrix}
$$

(https://books.physics.oregonstate.edu/GMM/kronecker.html). 

So to be more compact, they just write it as $I = \delta_{ij}$ sometimes. In physics, it becomes $\delta_i^j$, due to some common convention about what superscript and subscript indices mean.

Now go back to $\delta_{\alpha}^{\gamma} = \Lambda^{\bar{\beta}}_{\alpha} \Lambda^{\gamma}_{\bar{\beta}}$, we can see this is just saying $\Lambda^{\bar{\beta}}_{\alpha} \Lambda^{\gamma}_{\bar{\beta}}$ should be equal to an identity matrix.

Also note that the brackets in $\Lambda(\underset{\sim}{v})$ and $\Lambda(−\underset{\sim}{v})$ means they are parametrized by $\underset{\sim}{v}$ and $-\underset{\sim}{v}$, i.e. the brackets are not meant to emphasize a multiplication. And $\Lambda(\underset{\sim}{v})$ and $\Lambda(−\underset{\sim}{v})$ are inverse matrices.

## What does the eta (tensor) do in dot product

Related to this on the blackboard

$$
\vec{A} \cdot \vec{B} = (A^{\alpha} \vec{e_{\alpha}}) \cdot (B^{\beta} \vec{e_{\beta}}) \\
= A^{\alpha} B^{\beta} \;\; \vec{e_{\alpha}} \cdot \vec{e_{\beta}} \\
= A^{\alpha} B^{\beta} \;\; \eta_{\alpha \beta}
$$

**What alpha beta mean here**

They are **summation indices**, both running over the 4 components (0, 1, 2, 3) of the *same reference frame*.

That is:

* $\alpha$ indexes the component of $\vec A$,
* $\beta$ indexes the component of $\vec B$,
* both refer to the same set of basis vectors ${\vec e_\alpha}$.

So this is in one frame, one basis, but two dummy **summation** indices.

> The summation part is important. I was not used to thinking about it as summations yet at that time.

We could also pick any other 2 index labels and write

$$
\vec A \cdot \vec B = A^{\mu} B^{\nu} \; (\vec e_{\mu} \cdot \vec e_{\nu}),
$$

where the inner product of the basis vectors is what encodes the **metric**:

$$
\vec e_{\mu} \cdot \vec e_{\nu} = g_{\mu\nu}.
$$

Hence,

$$
\vec A \cdot \vec B = g_{\mu\nu} A^{\mu} B^{\nu}.
$$

**Why two indices appear**

It’s just because the product involves two vectors:

$$
(A^{\alpha} \vec e_\alpha) \cdot (B^{\beta} \vec e_\beta),
$$

so we need **one index for each expansion** - one for the A-expansion, one for the B-expansion.

When we actually carry out the dot product, we define

$$
\vec e_\alpha \cdot \vec e_\beta \doteq \eta_{\alpha\beta},
$$

and the result becomes

$$
A^{\alpha} B^{\beta} \eta_{\alpha\beta}.
$$

So we can think of the pair $(\alpha,\beta)$ as indexing the two slots of the metric tensor.

**Why the eta has those numbers**

$$
\eta_{\alpha\beta} = 
\begin{pmatrix}
-1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{pmatrix}
$$

Conceptually, $\eta_{\alpha\beta}$ is a tensor. But before the lecture gets to that, we can temp see it like a simple matrix here.

The bottom line is $\eta_{\alpha\beta}$ is the object that tells us how to take dot products in spacetime.

To break it down, first we need to know some rules for choosing coord system for spacetimes, in the realm of relativity.

In flat spacetime (special relativity), we usually choose a basis that is **orthonormal** with respect to the Minkowski metric - that is:

$$
\vec e_0 \cdot \vec e_0 = -1, \quad
\vec e_i \cdot \vec e_j = +\delta_{ij}, \quad
\vec e_0 \cdot \vec e_i = 0.
$$

Here $i,j = 1,2,3$ for the spatial directions. These look strage because the 1st one seems to be a self-dot-product that results in a negative number, and the 2nd and 3rd one also look like they have more meaning than pure math calculations.

This is all because linear algebra is all talking in the Euclidean space, but in relativity, spacetime is modeled as a pseudo-Euclidean (or Lorentzian) space. Formally:

$$
\text{signature of Minkowski metric}=(−,+,+,+).
$$

Therefore, for the time basis vector $e_0=(1,0,0,0)$ we pick:

$$e_0 \cdot e_0 = -1$$
. This just means we’re using a different notion of “distance.” (metric)

(https://archive.lib.msu.edu/crcmath/math/math/m/m264.htm)

> Remember the instructor's verbal tip on how topologies and geometries are separate things, and the donut vs coffee mug (with a handle) analogy in lecture 1?

Secondly, we notice that defining it like this works mathematically. Recall that we can write

$$
\vec A = (A^0, A^1, A^2, A^3), \\
\vec B = (B^0, B^1, B^2, B^3),
$$

. So then the $A^{\alpha} B^{\beta} \eta_{\alpha\beta}$ would be equivalent to $A^T \eta B$, and would give the distance in spacetime as we want:

$$
-A^0 B^0 + A^1 B^1 + A^2 B^2 + A^3 B^3
$$

. Just be ware that we need to get out of thinking about it in expanded matrix form ASAP, for our own benefits in the future. This here is just for easier understanding when we start our journey.

# Lecture 3

## Contravariant covariant “upstairs” and “downstairs”

Legacy terminology. For simplicity (but not take as truth) think about them like:

- Contravariant is like a "locational" vector field
- Covariant is like its gradient operators

, or:

- A contravariant vector's components change in the opposite way to the basis vectors, while
- a covariant vector's components change in the same way.

But in relativity, they are not very distinguishable. Quoting Schutz:

> These opposing transformation laws gave rise to the old names 'contravariant' and 'covariant'. What we call a vector was called contravariant because its components obey the law opposite ('contra') the basis vectors. Similarly one-forms were 'covariant vectors' because their components go with the basis vectors. The modern viewpoint emphasises the fact that neither the vector nor the one-form is in fact changed by a basis transformation: they are coordinate-independent geometrical objects. Therefore, modern terminology has dropped the old names because they emphasize the coordinate-dependent descriptions of these objects.

So in this course, the instructor uses “upstairs” and “downstairs”. In order to transition, we can think of their relation as

- Upstairs index: → contravariant component
- Downstairs index: → covariant component

In curved spacetime, coordinate transformations act differently on vectors and 1-forms. For a coordinate transformation ( x^\mu \to x'^\mu(x) ):

* Contravariant vector transforms as
  $$
  V'^{\mu} = \frac{\partial x'^{\mu}}{\partial x^{\nu}} V^{\nu}
  $$

* Covariant vector transforms as
  $$
  V'*{\mu} = \frac{\partial x^{\nu}}{\partial x'^{\mu}} V*{\nu}
  $$

Notice one uses the Jacobian, the other the *inverse* Jacobian.
That’s why they are different “types” of objects mathematically. But again, in relativity they are all seen as tensors.

## Simple contraction of tensors

### References

- Inner outer products: [phys.libretexts.org/.../19.06%3A_Appendix_-_Tensor_Algebra](https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/19%3A_Mathematical_Methods_for_Classical_Mechanics/19.06%3A_Appendix_-_Tensor_Algebra)
- Nice graphs: [https://www.ita.uni-heidelberg.de/.../tensor.pdf](https://www.ita.uni-heidelberg.de/~dullemond/lectures/tensor/tensor.pdf)
- Simple excercises: [https://grinfeld.org/books/An-Introduction-To-Tensor-Calculus/Chapter11.html](https://grinfeld.org/books/An-Introduction-To-Tensor-Calculus/Chapter11.html)
- Simple example with numbers: https://www.youtube.com/watch?v=OoT8kty3HPA

### Even simpler 2 dim tensors example

Let the contravariant components be:

$$
T^{\mu\nu} =
\begin{pmatrix}
T^{00} & T^{01} \\
T^{10} & T^{11}
\end{pmatrix}
=
\begin{pmatrix}
2 & -1 \\
3 & 4
\end{pmatrix}
$$

Let the covariant components (a different tensor) be:

$$
T_{\mu\nu} =
\begin{pmatrix}
5 & 2 \\
-3 & 1
\end{pmatrix}
$$

These two do **not** have any symmetry relationship.

Indices run over $\mu,\nu = 0, 1$.


By definition:

$$
T^{\mu\nu} T_{\mu\nu}
= T^{00}T_{00} + T^{01}T_{01} + T^{10}T_{10} + T^{11}T_{11}
$$

Now plug in the numbers:

1. $T^{00} T_{00} = 2 \cdot 5 = 10$
2. $T^{01} T_{01} = (-1) \cdot 2 = -2$
3. $T^{10} T_{10} = 3 \cdot (-3) = -9$
4. $T^{11} T_{11} = 4 \cdot 1 = 4$

Sum them:

$$
T^{\mu\nu}T_{\mu\nu}
= 10 - 2 - 9 + 4
= 3
$$

## Simple outer product example

The **outer product** takes two tensors and produces a new tensor with **more indices**.

Let

$$
v =
\begin{pmatrix}
2 \\
3
\end{pmatrix}, \quad
w =
\begin{pmatrix}
1 \\
4
\end{pmatrix}
$$

Then the outer product $v \otimes w$ is:

$$
(v \otimes w)^{\mu\nu}
=

\begin{pmatrix}
2 \cdot 1 & 2 \cdot 4 \\
3 \cdot 1 & 3 \cdot 4
\end{pmatrix}
=

\begin{pmatrix}
2 & 8 \\
3 & 12
\end{pmatrix}
$$

It is just multiplying independently across dimensions.

This is **not** a matrix multiplication - it’s building a tensor by combining components of two vectors.

If you have:

* $A^{\mu_1 \cdots \mu_r}$: a rank-r tensor
* $B^{\nu_1 \cdots \nu_s}$: a rank-s tensor

Their outer product is a rank $r+s$ tensor:

$$
(A \otimes B)^{\mu_1 \cdots \mu_r \nu_1 \cdots \nu_s}
=

A^{\mu_1 \cdots \mu_r} , B^{\nu_1 \cdots \nu_s}
$$

## An index lowering example

In relativity, we often use the **Minkowski metric** for flat spacetime. It's a fundamental tool for calculations in Special Relativity. We'll use this metric in our example.

### The Metric Tensor ($g_{\mu\nu}$)

The metric tensor defines the geometry of spacetime. For the flat spacetime of Special Relativity, we use the Minkowski metric, commonly denoted as $\eta_{\mu\nu}$ or $g_{\mu\nu}$. We'll use the convention with the signature `(-, +, +, +)`, which is what is used in lecture 1-3.

> Signature (+, -, -, -) also seems to be commonly used.

The indices $\mu$ and $\nu$ run from 0 to 3, where 0 is the time component and 1, 2, and 3 are the spatial components.

### The Index Lowering Operation

We want to lower the second index, $\nu$, to obtain the mixed tensor $T^\mu_\nu$. The rule for lowering an index is to contract the tensor with the metric tensor. The formula is:

$T^\mu_\nu = g_{\nu\sigma} T^{\mu\sigma}$

This formula uses the **Einstein summation convention**, which means that whenever an index appears once as an upper index and once as a lower index in a term, you should sum over all possible values of that index (from 0 to 3). In this case, the index $\sigma$ is the "dummy" summation index.

Expanded, the formula for a single component of the new tensor looks like this:

$T^\mu_\nu = \sum_{\sigma=0}^{3} g_{\nu\sigma} T^{\mu\sigma} = g_{\nu0}T^{\mu0} + g_{\nu1}T^{\mu1} + g_{\nu2}T^{\mu2} + g_{\nu3}T^{\mu3}$

In terms of matrix multiplication, this operation is equivalent to multiplying the matrix for $T^{\mu\nu}$ by the matrix for $g_{\mu\nu}$. However, be careful with the order and the indices. The formula $T^\mu_\nu = g_{\nu\sigma} T^{\mu\sigma}$ shows that for each row $\mu$ of $T$, the new row is a linear combination of the old row's components, with the coefficients given by the row of the metric that corresponds to the index $\nu$.

A simpler way to see this is that we are multiplying the columns of $T^{\mu\nu}$ by the corresponding diagonal elements of $g_{\mu\nu}$.

### Calculation of it

*   **Metric Tensor ($g_{\mu\nu}$):**
    $g_{\mu\nu} = \begin{pmatrix} -1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$

*   **Upstairs (contravariant) Tensor ($T^{\mu\nu}$):**
    $T^{\mu\nu} = \begin{pmatrix} 2 & 0 & 1 & -1 \\ -1 & 0 & 3 & 2 \\ -1 & 1 & 0 & 0 \\ -2 & 1 & 1 & -2 \end{pmatrix}$

*   **Index Lowering ($T^\mu_\nu = g_{\nu\sigma} T^{\mu\sigma}$):** This operation multiplies the columns of $T^{\mu\nu}$ by the corresponding diagonal elements of our new $g_{\mu\nu}$.

    *   The first column ($\nu=0$) gets multiplied by **-1**.
    *   The second column ($\nu=1$) gets multiplied by **+1**.
    *   The third column ($\nu=2$) gets multiplied by **+1**.
    *   The fourth column ($\nu=3$) gets multiplied by **+1**.

*   **The Resulting Mixed Tensor ($T^\mu_\nu$):**

    $T^\mu_\nu = \begin{pmatrix}
    -1 \cdot (2) & 1 \cdot (0) & 1 \cdot (1) & 1 \cdot (-1) \\
    -1 \cdot (-1) & 1 \cdot (0) & 1 \cdot (3) & 1 \cdot (2) \\
    -1 \cdot (-1) & 1 \cdot (1) & 1 \cdot (0) & 1 \cdot (0) \\
    -1 \cdot (-2) & 1 \cdot (1) & 1 \cdot (1) & 1 \cdot (-2)
    \end{pmatrix}$

    So, with the `(-, +, +, +)` signature, the final mixed tensor is:

    $T^\mu_\nu = \begin{pmatrix} -2 & 0 & 1 & -1 \\ 1 & 0 & 3 & 2 \\ 1 & 1 & 0 & 0 \\ 2 & 1 & 1 & -2 \end{pmatrix}$

The process of "raising" an index is similar but uses the inverse metric, $g^{\mu\nu}$. For the Minkowski metric, the inverse metric happens to be identical to the metric itself.

# Lecture 4

## “Trajectory” of photon is null curve

### Intuition and proper time

In general relativity, the motion of any particle is described by its **4-momentum** or **4-velocity**. For a massive particle, its 4-velocity satisfies

$$
u^\mu u_\mu = -1
$$

(with the (-,+,+,+) signature) - this means its worldline is **timelike**.

Let's look at the 4-velocity in differential form. For a "massive" particle it is:

$$
u^\mu = \frac{dx^\mu}{d\tau},
$$

where $\tau$ is **proper time**, the time measured by a clock moving with the particle. This is the *only* time parameter that makes sense for defining a velocity in spacetime.

Proper time satisfies:

$$
d\tau^2 = - \frac{1}{c^2} ds^2.
$$

For a timelike worldline (massive particle):

$$
ds^2 < 0 \quad\Rightarrow\quad d\tau^2 > 0.
$$

For a null worldline (photon):

$$
ds^2 = 0 \quad\Rightarrow\quad d\tau = 0.
$$

This is the key: A photon experiences no proper time.

**Because photons travel at the speed of light.**

Reminder about timelike and such:

| Type of vector       | Condition         | Physical meaning                           |
| -------------------- | ----------------- | ------------------------------------------ |
| **Timelike**         | $v^\mu v_\mu < 0$ | Massive particles (slower than light)      |
| **Null (lightlike)** | $v^\mu v_\mu = 0$ | Light or anything traveling at exactly (c) |
| **Spacelike**        | $v^\mu v_\mu > 0$ | Cannot be a physical trajectory            |

So anyway:

$$
u^\mu = \frac{dx^\mu}{d\tau}
$$

is undefined for photons because the denominator is zero.

### The physical reasoning

In special relativity, the invariant interval is

$$
ds^2 = -c^2 dt^2 + dx^2 + dy^2 + dz^2.
$$

For a light ray:

$$
ds^2 = 0.
$$

This means the spacetime separation between any two points on a photon's trajectory is zero.

In GR, this generalizes to:

$$
g_{\mu \nu} ,\frac{dx^\mu}{d\lambda} ,\frac{dx^\nu}{d\lambda} = 0,
$$

where $\lambda$ is any affine parameter along the photon’s worldline.

Thus the **tangent vector** to a photon’s trajectory has zero norm.

Since 4-momentum is proportional to the tangent vector for a geodesic:

$$
p^\mu p_\mu = 0.
$$

## Normalizing 4-velo?

Normalizing the 4-velocity uses the fact that $d\tau$ is the “clock tick”

For massive particles:

$$
u^\mu u_\mu = g_{\mu\nu} \frac{dx^\mu}{d\tau} \frac{dx^\nu}{d\tau}
= \frac{ds^2}{d\tau^2} = -c^2.
$$

So conventionally we divide by $c$ and write:

$$
u^\mu u_\mu = -1.
$$

For photons, $d\tau = 0 → u^\mu$ is not defined.

## How mass connects to the norm

Start from the SR energy-momentum relation:

$$
E^2 = (pc)^2 + (mc^2)^2.
$$

This is the “expanded” form. There is a covariant form:

$$
p^\mu p_\mu = - m^2 c^2,
$$

We can confirm the covariant form is equivalent to the expanded form. Start by expressing it in 4-vector language. Define the **4-momentum**:

$$
p^\mu = \left(\frac{E}{c}, p_x, p_y, p_z\right).
$$

The Minkowski metric (in −,+,+,+ signature):

$$
\eta_{\mu\nu} = \text{diag}(-1,1,1,1).
$$

Compute the invariant:

$$
p^\mu p_\mu
= \eta_{\mu\nu} p^\mu p^\nu
= -\left(\frac{E}{c}\right)^2 + p^2.
$$

So:

$$
p^\mu p_\mu = -\left(\frac{E}{c}\right)^2 + p^2.
$$

Set this equal to the invariant $-m^2 c^2$:

$$
-\left(\frac{E}{c}\right)^2 + p^2 = -m^2 c^2.
$$
Now multiply both sides by $(-c^2)$:

$$
E^2 - p^2 c^2 = m^2 c^4.
$$

Rearrange:

$$
E^2 = p^2 c^2 + m^2 c^4.
$$

So they are the same relation. Then from the covariant form we see:

* If $m > 0$: $p^\mu p_\mu < 0$ (timelike)
* If $m = 0$: $p^\mu p_\mu = 0$ (null)
* If $m^2 < 0$: impossible

This is why the norm of the 4-momentum - and by extension, the norm of the 4-velocity - tells us the rest mass.

Indeed:

$$
u^\mu u_\mu = -1 \quad\Leftrightarrow\quad m \neq 0.
$$

Another example is a photon has **0 rest mass**.

## What is h omega

The instructor wrote this for a photon:

$$
p^\mu = \bar{h}\omega (1, \hat{1})
$$

where:

* $p^\mu$ is the photon **4-momentum**,
* $\hat{1}$ is a **unit spatial direction vector** (a unit 3-vector),
* $\omega$ is the **angular frequency of the photon**,
* $\bar{h}$ is the **reduced Planck constant** (better symbol is $\hbar$).


$\hbar \omega$ is literally the **energy** of a photon.

$$
E = \hbar \omega.
$$

This is the quantum-energy formula relating:

* $\omega$ (angular frequency)
* $E$ (energy)
* $\hbar$ (reduced Planck constant)

Angular frequency $\omega$ is related to wavelength and ordinary frequency:

$$
\omega = 2\pi f = \frac{2\pi c}{\lambda}.
$$

So:

$$
E = \hbar \omega = h f.
$$

Why write this with $\hbar \omega$ as an example?

Maybe to quickly demonstrate 4-momentum of a massless particle. For a massless particle the 4-momentum must satisfy:

$$
p^\mu p_\mu = 0.
$$

And

$$
p^\mu = (E/c, \vec{p}).
$$

A photon has:

$$
|\vec{p}| = \frac{E}{c}.
$$

So the simplest way to write a null 4-vector is:

$$
p^\mu = E (1, \hat{n}),
$$

or if using $\hbar\omega$ instead of E:

$$
p^\mu = \hbar \omega (1, \hat{n}).
$$

This works because:

$$
p^\mu p_\mu = -(\hbar\omega)^2 + (\hbar\omega)^2 = 0.
$$

So the expression

$$
  p^\mu = \hbar \omega (1, \hat{1})
$$

is a compact way to express the **null 4-momentum** of a photon.

## Levi–Civita and Permutation Parity

For (n)-dimensional space (most commonly (n=3)):

$$
\varepsilon_{i_1 i_2 i_3} =
\begin{cases}
+1 & \text{if }(i_1,i_2,i_3)\text{ is an even permutation of }(1,2,3), \\
-1 & \text{if }(i_1,i_2,i_3)\text{ is an odd permutation of }(1,2,3), \\
0 & \text{if any indices are repeated.}
\end{cases}
$$

So the Levi–Civita symbol is essentially a “parity detector.”

What Are Even and Odd Permutations Here?

Take the ordered tuple $(1,2,3)$.

A permutation like:

* $(1,2,3)$ itself is an **even** permutation → $ \varepsilon_{123} = +1 $
* $(2,3,1)$ (cyclic shift) is also **even** → $ \varepsilon_{231} = +1 $
* $(3,1,2)$ is **even** → $ \varepsilon_{312} = +1 $

These three are even because each can be made from **2 transpositions** (or 0 in the trivial case).

Examples of odd permutations:

* $(1,3,2)$ → $ \varepsilon_{132} = -1 $
* $(3,2,1)$ → $ \varepsilon_{321} = -1 $
* $(2,1,3)$ → $ \varepsilon_{213} = -1 $

These require an odd number of swaps to get from $(1,2,3)$.

**Quick Example: Cross Product**

$$
(a \times b)_i = \sum_{jk} \varepsilon_{ijk} a_j b_k
$$

If $(i,j,k)$ is an even permutation of $(1,2,3)$, the component matches the usual right-hand rule.
If odd, it flips sign.

**Geometric Meaning**

* Even permutations preserve orientation
  → same handedness → (+1)

* Odd permutations reverse orientation
  → flip handedness → (-1)

If we swap two basis vectors in a right-handed coordinate system, we get a left-handed one, and the Levi–Civita symbol reflects this sign flip.

## Derive 4D Gauss Thm thru a V4 integral

Related to this on the black board.

...

$$
\begin{aligned}
\int_{\partial V^4} N^{\alpha} d\Sigma_{\alpha}
&=
\underbrace{\int_{t=t_2} N^0 dxdydz}_{\text{upper time face}}
-\underbrace{\int_{t=t_1} N^0 dxdydz}_{\text{lower time face}}\\
&\quad +\underbrace{\int_{x=x_2} N^1 dtdydz}_{\text{face at }x_2}
-\underbrace{\int_{x=x_1} N^1 dtdydz}_{\text{face at }x_1}\\
&\quad +\underbrace{\int_{y=y_2} N^2dtdxdz}_{\text{face at }y_2}
-\underbrace{\int_{y=y_1} N^2 dtdxdz}_{\text{face at }y_1}\\
&\quad +\underbrace{\int_{z=z_2} N^3 dtdxdy}_{\text{face at }z_2}
-\underbrace{\int_{z=z_1} N^3 dtdxdy}_{\text{face at }z_1}
=0 .
\end{aligned}
$$

Let $t_2 \to t_1​+dt$, rearrange:

$$
\begin{aligned}
\int_{t_1+dt}N^0dxdydz-\int_{t_1}N^0dxdydz
= -dt
\Big[
\int_{x_2}N^1dydz-\int_{x_1}N^1dydz \\
+\int_{y_2}N^2dxdz-\int_{y_1}N^2dxdz \\
+\int_{z_2}N^3dxdy-\int_{z_1}N^3dxdy
\Big].
\end{aligned}
$$

divide by $dt$, next going to take limit

...

**Why can $dt$ be moved out of spatial integrals**

Over an infinitesimally thin time slab $t \in [t_1, t_1+dt]$, the integral over $t$ behaves like:

$$
\int_{t_1}^{t_1+dt} N^1(t,y,z), dt
= N^1(t_1,y,z), dt + O(dt^2)
$$

This uses the fact that *any smooth function can be approximated by its value at one point* over an infinitesimal interval, we get:

$$
N^1(t,y,z) = N^1(t_1,y,z) + (t - t_1)\frac{\partial N^1}{\partial t}(t_1,y,z) + \cdots
$$

When integrating from $t_1$ to $t_1+dt$:

* The first term gives $N^1(t_1,y,z), dt$.
* The next term gives $\frac{1}{2}\partial_t N^1(t_1,y,z) (dt)^2$ $\to O((dt)^2)$.
* Higher terms give even smaller orders.

Next going back to the hypersurface measure. For the face $x=\text{const}$ it is:

$$
d\Sigma_1 = dtdydz.
$$

Use it on the right side and the equation becomes:

$$
\int_{x=x_2} N^1 dt dy dz
= \int_{y_1}^{y_2} \int_{z_1}^{z_2}
\left(\int_{t_1}^{t_1+dt} N^1(t,y,z) dt \right) dydz.
$$

Now expand the inner integral using that approximation we got, we obtain:

$$
\int_{x=x_2} N^1  dt dy dz
= dt \int_{y_1}^{y_2} \int_{z_1}^{z_2} N^1(t_1,y,z) dy dz + O(dt^2).
$$

Now factoring out the $dt$ is possble. Same for other spatial integrals.