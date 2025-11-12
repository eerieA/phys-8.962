<!-- TOC -->

- [Lecture 1](#lecture-1)
    - [Deriving the Lorentz Transformations](#deriving-the-lorentz-transformations)
- [Lecture 2](#lecture-2)
    - [Kronecker delta and how to use for basis vecs](#kronecker-delta-and-how-to-use-for-basis-vecs)
    - [What does the eta (tensor) do in dot product](#what-does-the-eta-tensor-do-in-dot-product)

<!-- /TOC -->

# Lecture 1

## Deriving the Lorentz Transformations

Tutorial on deriving it for $v$ along $x$ axis:  
https://www.wetube.com/watch?v=FvqutkaPmas

> I like the explaination on why the Lorentz Transformation was chosen to be linear. This instructor thought it is because of the principle of relativity:
>
> - The transformation should be "reversible".
> - Also the transformation from $K$ to $K'$ and its "reverse" should be of the same form (no "previledged" IRF needed to determine which is which), and a linear one satisfies that.

Tutorial on deriving it for $v$ along an arbitrary axis:  
https://www.wetube.com/watch?v=Afd34FuG65A


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

In flat spacetime (special relativity), we usually choose a basis that is **orthonormal** with respect to the Minkowski metric — that is:
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