# Nelder-Mead-Downhill-Simplex

## Introduction
The Nelder–Mead simplex method [49, 50] is an algorithm for finding the local minima of a function $f:R^n→R, called the objective function in optimization parlance. The idea
behind it is very geometric. It starts with a suitably chosen initial simplex64 in Rn and
evolves it according to certain rules. The simplex undergoes iterative expansions, contrac-
tions, and displacements, and if everything goes according to plan, it moves “downhill”
toward a local minimum of f . The term “downhill simplex” was used in [53] to refer to
this behavior. A significant merit of the algorithm is that it is gradient-free, that is, it does
not require a knowledge of the function’s gradient, which may be difficult or expensive
to compute. Consequently, it may be applied to minimize nondifferentiable functions,
but see the caveats in Section 18.3.
18.2 The algorithm
The algorithm works in iterative cycles. I will explain what is done in one cycle. The
cycles are repeated until certain convergence criteria are met. The calculations involve
four parameters, β r , βe , βc , and β s , which permit some customization of the algorithm,
although almost always these are taken as
β r = 1,
βe = 2,
1
βc = ,
2
1
βs = .
2
Now, without further ado, let us begin with the description of the algorithm.
Rank the vertices. Let us write x(i ) ∈ Rn , i = 0, 1, . . . , n, for the simplex’s n + 1 vertices
at the start of a cycle, and let yi = f (x(i ) ) be the function values at those vertices.
Three of the n + 1 vertices play key roles in what happens next. The best vertex is
64 A simplex in n dimensions is the generalization of a triangle (two dimensions) and a tetrahedron (three
dimensions). It is the convex hull of n + 1 points in the n-dimensional Euclidean space. Equivalently, it
may be thought of as the “solid” object formed by a set of n vectors emanating from a common point in the
n-dimensional space. If the vectors form a linearly independent set, then the simplex has a positive n-dimensional
volume; otherwise it is a degenerate simplex.
193194
Chapter 18. The Nelder–Mead downhill simplex
one with the lowest function value (it’s best because we are aiming to minimize the
function). The worst vertex is one with the highest function value. The next to worst
vertex is one with the highest function value after excluding the worst vertex. We
write ia, i y, and i z for the indices of the best, next to worst, and worst vertices,
respectively, because a, y, and z are the first, next to last, and last letters of the
alphabet.65
The definitions of ia, i y, and i z are not quite unambiguous. What are the best,
next to worst, and worst vertices if the vector of the function values is [7, 7, 7, 7] or
[5, 5, 1, 1, 7, 7]? We partially remove the ambiguity by requiring that ia, i y, and i z
be distinct. Stated formally, we require that
ia = i y,
yi a = min yi ,
i
i y = i z,
i z = ia,
yi z = max yi ,
yi y = max yi
if n ≥ 2.
(18.1a)
i =i z
i
Thus, if the function values are [7, 7, 7, 7], then ia = 0, i z = 1, i y = 2 will do,
but ia = i z = i y = 0 is not acceptable. Similarly, if the function values are
[5, 5, 1, 1, 7, 7], then ia will have to be one of 2 or 3; and either of the choices i z = 4,
i y = 5 or i z = 5, i y = 4 will do.66
In defining the best, next to worst, and worst vertices above I have assumed implic-
itly that the simplex has at least three vertices, that is, the dimension of space n ≥ 2.
The one-dimensional case, n = 1, is special; the simplex is a line segment, and hence
it has only two vertices. From the definition of i y it follows that i y = ia. Therefore
picking distinct ia, i y, and i z is impossible. We make an exception in the case of
n = 1 and change the first line of the conditions (18.1a) to
ia = i y = i z
if n = 1.
(18.1b)
The ambiguities noted above are not the only ones. There are quite a few decisions
in the rest of the algorithm that call for tie-breakings of various sorts, e.g., making
comparisons with “<” versus “≤”. I will gloss over the precise tie-breaking rules—as
did Nelder and Mead in their paper—therefore your program may not work exactly
as someone else’s that breaks the ties differently. See [40] for a set of precise tie-
breaking rules that lead to a uniquely defined Nelder–Mead algorithm.
Reflect. Compute the centroid x̂ of the simplex face opposite to x(i z) ,
x̂ =
1
n i =i z
x(i ) ,
(18.2)
compute the “reflection” x(r ) of x(i z) through x̂,
x(r ) − x̂ = −β r (x(i z) − x̂),
(18.3a)
and let y (r ) = f (x(r ) ). If β r = 1, as it almost always is, then x(r ) is the true reflection
of x(i z) through x̂ because then x(r ) − x̂ = x̂ − x(i z) . The reflection diagram in
Figure 18.1 illustrates this.
At this point the algorithm branches into four cases, depending on the relationship
between the value of y (r ) and the vertex values yi , i = 0, 1, . . . , n.
65 Here I am breaking away from the tradition in mathematics of using single letters for variable names because
I could not think of suitable replacements for the double-letter names i a, i y, and i z. From the programming
point of view, however, these are quite natural; programmers routinely use multiletter symbols for variables.
66
Here, and everywhere else in this book for that matter, array indices begin with zero.18.2. The algorithm
195
x(i z)
x(i z)
x̂
x̂
reflection
x(r )
x(i z)
expansion
x(r )
x(i z)
x(e)
x(c)
x̂
x(c)
x̂
inner contraction
outer contraction
x(i a)
shrink
x(r )
Figure 18.1: These diagrams show all possible transformations of a Nelder–Mead sim-
plex in two dimensions (n = 2) with the standard choices of the parameter
values: β r = 1, βe = 2, βc = 1/2, β s = 1/2.
Case 1: y (r ) < yi a
This is a terrific occurrence; what used to be the worst vertex changes to the
best vertex upon reflection. Thus, the line that connects x(i z) through x̂ to
x(r ) appears to be a direction of steep descent of the function f . Take advan-
tage of the opportunity, and expand the simplex in that direction by a factor
of βe . Thus, calculate a point x(e) such that
x(e) − x̂ = βe (x(r ) − x̂),
(18.3b)
and let y (e) = f (x(e) ). Usually one takes βe = 2. This is illustrated in the
expansion diagram in Figure 18.1. Replace the vertex x(i z) with the better of
x(e) and x(r ) , and end the iteration cycle.
Remark 18.1. This is slightly different from Nelder and Mead’s original
prescription which is the following: Pick x(e) if y (e) < yi a ; else pick x(r ) .
Case 2: yi a ≤ y (r ) < yi y
Replace the vertex x(i z) with x(r ) , and end the iteration cycle.
Case 3: yi y ≤ y (r ) < yi z
Replacing x(i z) by x(r ) under these conditions will constitute some improve-
ment, but it will not be very useful since the worst vertex remains the worst
vertex after the replacement. Thus, try another tack: push x(r ) toward the196
Chapter 18. The Nelder–Mead downhill simplex
centroid x̂ by a factor of βc to arrive at a point x(c) ,
x(c) − x̂ = βc (x(r ) − x̂),
(18.3c)
and let y (c) = f (x(c) ). Usually one takes βc = 1/2. The outer contraction
diagram in Figure 18.1 illustrates this. If y (c) < y (r ) , accept x(c) as a replace-
ment for the vertex x(i z) and end the cycle; otherwise shrink the simplex (see
below) and end the cycle.
Case 4: yi z ≤ y (r )
In this case a reflection is of no use at all since the worst point becomes even
worse. Try another tack: push x(i z) toward the centroid x̂ by a factor of βc
to arrive at a point x(c) ,
x(c) − x̂ = βc (x(i z) − x̂),
(18.3d)
and let y (c) = f (x(c) ). The inner contraction diagram in Figure 18.1 illustrates
this. If y (c) < yi z , accept x(c) as a replacement for the vertex x(i z) and end the
cycle; otherwise shrink the simplex (see below) and end the cycle.
Shrink. The shrink process referred to in cases 3 and 4 above shrinks the simplex toward
its best vertex x(i a) by moving all other vertices toward x(i a) by a factor of β s :
)
x(inew
− x(i a) = β s (x(i ) − x(i a) ) for all i = ia.
(18.3e)
Usually one takes β s = 1/2. This is illustrated in the shrink diagram in Figure 18.1.
Are we there yet? There is no universally accepted stopping criterion for the Nelder–
Mead algorithm. Nelder and Mead’s original article [49] suggested the criterion
'
( 
(1 n
)
(y − ȳ)2 < ε,
n i =0 i
motivated by applications in statistics, where ȳ is the average of the yi ’s, and ε is a
prescribed measure of accuracy. The stopping criterion proposed in [53] is
2
|yi z − yi a |
|yi z | + |yi a |
<ε
for some prescribed ε. Either criterion is met if the variation of the objective func-
tion’s values over the simplex’s vertices is small. The simplex’s size does not enter
into the consideration.
Wright [80] quotes others who have suggested
max x(i ) − x(i a) ≤ ε max(1, x(i a) )
i =i a
for the stopping criterion, where ε is a prescribed measure of accuracy. This crite-
rion is met if the simplex is sufficiently small. The objective function’s values do
not enter into the consideration.
The following hybrid criterion accounts for both the simplex’s size and the objec-
tive function’s values. It is quite robust and has worked well for my purposes, so
that’s what I prescribe for our program.18.3. Problems with the Nelder–Mead algorithm
197
The user supplies a number, h, that establishes a “length scale” in Rn . A good h
should be of an order of magnitude commensurate with the features of interest in
the objective function’s domain. The length scale h serves dual purposes. For one
thing, it enters the stopping criterion described below, and for another thing, it is
used to set the simplex’s initial size.
The user also supplies a dimensionless number, τ (for tolerance), that determines the
algorithm’s accuracy. The program lets ε = h × τ and declares that the algorithm
has converged if
x(i z) − x(i a) ≤ ε and
|yi z − yi a | ≤ ε2 .
(18.4)
This assumes that ε is suitable for measuring the “smallness” of both x and y quan-
tities. It also assumes that the objective function behaves like a quadratic at the
minimizing point, which is generally the case if the function is smooth. You are
free to replace these generic conditions with others which you deem better suited
to problems of special interest to you.
To avoid runaway computations, should the algorithm fail to converge, the program
keeps a count of the number of times it evaluates the objective function. It halts the
iterations if the evaluation count exceeds a user-specified ceiling.