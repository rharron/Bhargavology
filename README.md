# Bhargavology
Sage (python) scripts for working with cubic and quartic rings via binary cubic forms and pairs of ternary quadratic forms in the spirit of Delone–Faddeev, Davenport–Heilbronn, Gan–Gross–Savin, Bhargava, et al.

**Note:** this code is still in development. I'm releasing it because I think it can be useful to people out there despite this fact.

### Table of contents

* [Basic usage comments](#basic-usage-comments)
* [Cubic rings and binary cubic forms](#cubic-rings-and-binary-cubic-forms)
  
  * [Some basic functionality](#some-basic-functionality)
  * [Dedekind's example](#dedekinds-example)
  * [Other functions of interest](#other-functions-of-interest)

## Basic usage comments

To start using this code in Sage use:
```
sage: load('PairOfTernaryQuadraticForms.sage')
```

This loads the scripts for pairs of ternary quadratic forms, ternary quadratic forms, binary cubic forms, and binary quadratic forms. If you won't be dealing with quartic rings or pairs of ternary quadratic forms, you can simply do:

```
sage: load('BinaryCubicForm.sage')
```

This will load the bianary cubic form and binary quadratic form code.

## Cubic rings and binary cubic forms

In this section, I give some examples of how to use my code. I've tried to include some interesting examples when possible. For mathematical details on the relation between cubic rings and binary cubic forms, see Sections 2 & 3 of [Bhargava–Shankar–Tsimerman](http://dx.doi.org/10.1007/s00222-012-0433-0).

### Some basic functionality

The simplest way to create a binary cubic form is to hand it a list of 4 integers:

```
sage: BCF = BinaryCubicForm([1, -1, -2, -8]); BCF
X^3 - X^2*Y - 2*X*Y^2 - 8*Y^3
```

But we can also work over other base rings, and we can customize the variable names:

```
sage: K.<z> = GF(4)
sage: BCF = BinaryCubicForm([z, -1, -1, z+1], base_ring=K, var_names='x, y'); BCF
(z)*x^3 + x^2*y + x*y^2 + (z + 1)*y^3
```

The code can compute discriminants of and actions on binary cubic forms:

```
sage: R.<a, b, c, d> = PolynomialRing(QQ)
sage: BCF = BinaryCubicForm([a, b, c, d], base_ring=R)
sage: BCF.disc()
b^2*c^2 - 4*a*c^3 - 4*b^3*d + 18*a*b*c*d - 27*a^2*d^2
sage: R2.<N> = PolynomialRing(R)
sage: BCF.action(Matrix(2, [1,N,0,1]))
(d*N^3 + c*N^2 + b*N + a)*X^3 + (3*d*N^2 + 2*c*N + b)*X^2*Y + (3*d*N + c)*X*Y^2 + d*Y^3
```

Note that the action is via the so-called twisted action:
```
sage: BinaryCubicForm([1, -1, -2, -8]).action(Matrix(2, [1,0,0,1/2]))
2*X^3 - X^2*Y - X*Y^2 - 2*Y^3
```

You can ask for a multiplication table for the normalized basis:
```
sage: R.<a, b, c, d> = PolynomialRing(QQ)
sage: BCF = BinaryCubicForm([a, b, c, d], base_ring=R)
        | 1       omega                        theta
+-------+-------+----------------------------+----------------------------+
  1     | 1       omega                        theta
  omega | omega   (-b)*omega + a*theta - a*c   -a*d
  theta | theta   -a*d                         (-d)*omega + c*theta - b*d
sage: BinaryCubicForm([0,-1,1,0]).multiplication_table()
        | 1       omega   theta
+-------+-------+-------+-------+
  1     | 1       omega   theta
  omega | omega   omega   0
  theta | theta   0       theta
```
This latter multiplication table is that of the cubic ring **Z**<sup>3</sup> with basis
((1,1,1), (0,1,0), (0,0,1)). We can't quite create the corresponding cubic ring in Sage,
but the Bhargavolgy code does the next best thing, it creates a cubic **Q**-algebra with
that basis:

```
sage: BCF0 = BinaryCubicForm([0,-1,1,0])
sage: A = BCF0.algebra();
sage: one, omega, theta = A.basis()
sage: (omega+theta).characteristic_polynomial()
x^3 - 2*x^2 + x
```

When a binary cubic form corresponds to an order in a number field (i.e. when the form
is irreducible), we can get the associated order:
```
sage: BCF = BinaryCubicForm([1, -1, -2, -8])
sage: BCF.order()
Order in Number Field in eta with defining polynomial x^3 - x^2 - 2*x - 8
```

For more examples of usage, see the examples below, as well as the examples in the code itself.

### Dedekind's example

One major advantage of using binary cubic forms is that it allows you to work with non-monogenic cubic rings. In this example, I'll use Dedekind's original example of a non-monogenic ring of integers to illustrate some functionality of this Bhargavology code. Dedekind's example is **Q**(&alpha;) where &alpha; is a root of *f*(*x*)&nbsp;=&nbsp;*x*<sup>3</sup>&minus;*x*<sup>2</sup>&minus;2*x*&minus;8.

Let's get the order corresponding to *f*(*x*) and see that it's not maximal:

```
sage: x = polygen(ZZ)
sage: BCF = BinaryCubicForm([1, -1, -2, -8]); BCF
X^3 - X^2*Y - 2*X*Y^2 - 8*Y^3
sage: O = BCF.order()
sage: O.is_maximal()
False
```

Since BCF(*x*, 1) = *f*(*x*), *O*=**Z**[&alpha;]. We can ask the Bhargavology code for a binary cubic form corresponding to the maximal order containing *O* as follows:

```
sage: BCFmax = BCF.find_maximal_form(); BCFmax
2*X^3 - X^2*Y - X*Y^2 - 2*Y^3
sage: Omax = BCFmax.order()
sage: Omax.is_maximal()
True
```
We could have also gotten this "by hand" as follows. Factoring the discriminant of BCF, we can see that the only square dividing it is 2<sup>2</sup>, so if *O* is not maximal, then it has index 2 in the maximal order. Because the last coefficient of BCF is divisible by 2<sup>2</sup> and the second-to-last by 2, finding a form corresponding to a ring containing *O* with index 2 is straightforward:

```
sage: BCF.disc().factor()
-1 * 2^2 * 503
sage: BCF.action(Matrix(2, [1,0,0,1/2]))
2*X^3 - X^2*Y - X*Y^2 - 2*Y^3
```

We can see the multiplication table of *O* in the normalized basis 1, &omega;, &theta; (see section 2 of Bhargava–Shankar–Tsimerman, where they use the term "normal basis"):
```
sage: BCF.multiplication_table()
        | 1       omega               theta
+-------+-------+-------------------+-----------------------+
  1     | 1       omega               theta
  omega | omega   omega + theta + 2   8
  theta | theta   8                   8*omega - 2*theta - 8
```

Now that we have a binary cubic form corresponding to the maximal order, we can use it to understand some algebraic number theory. The cubic ring corresponding to a binary cubic form is monogenic if and only if the form represents &pm;1. We can thus see that the ring of integers of **Q**(&alpha;) is not monogenic by, say, reducing mod 2 and seeing that BCFmax doesn't represent 1 mod 2:

```
sage: [BCFmax(i, j) for i in GF(2) for j in GF(2)]
[0, 0, 0, 0]
```

The rest of this example doesn't really use the Bhargavology code, but I wanted to give an idea of what you can do with binary cubic forms. Dedekind's example is one where 2 always divides the index of any monogenic order, but using binary cubic forms, you can still determine the prime factorization of 2 in a manner analogous to Dedekind's theorem for polynomials (see Section 3 of del Corso–Dvornicich–Simon's "Decomposition of primes in non-maximal orders" for the mathematical details). Specifically, factoring BCFmax mod 2 gives polynomials that you should use to generate the prime ideals dividing 2*O*<sub>max</sub>:

```
sage: K = O.fraction_field()
sage: a, b, c, d = BCFmax.coefficients()
sage: alpha = BCFmax(x, 1).roots(K)[0][0]
sage: B = K.ideal([a, a*alpha+b, a*alpha^2+b*alpha+c]) #An ideal used for all factorizations
sage: facts = BCFmax.change_ring(GF(2)).polynomial().factor(); facts
Y * X * (X + Y)
sage: Ps = []
sage: for F, e in facts:
          FZ = F.change_ring(ZZ)
          Ps.append(K.ideal(2) + B^(e * F.degree()) * FZ(alpha, 1))
sage: map(lambda x : x.is_prime(), Ps)
[True, True, True]
sage: prod(Ps)
Fractional ideal (2)
```

The last two commands verify that we generated 3 prime ideals whose product is 2*O*<sub>max</sub>.

### Other functions of interest ###
In this section, I briefly mention the some other functions in the binary cubic forms code. For the full list, see the source code itself.

You can compute the Hessian (the trace zero form is 6 times the Hessian)
```
sage: BCF = BinaryCubicForm([0, -1, 1, 0])
sage: BCF1 = BinaryCubicForm([1, 0, 0, -10])
sage: R.<a, b, c, d> = PolynomialRing(QQ)
sage: BCF2 = BinaryCubicForm([a, b, c, d], base_ring=R)
sage: BCF.hessian()
X^2 - X*Y + Y^2
sage: BCF1.hessian()
90*X*Y
sage: BCF2.hessian()
(b^2 - 3*a*c)*X^2 + (b*c - 9*a*d)*X*Y + (c^2 - 3*b*d)*Y^2
sage: BCF1.trace_zero_form()
540*x*y
```

You can also compute the shape (which is a binary quadratic form over
different rings depending on what binary cubic form you have):
```
sage: BCF.shape()
X^2 - X*Y + Y^2
sage: BCF1.shape()
125.3228985075451?*x^2 + 581.697366308609?*y^2
sage: BCF2.shape(check=False, use_AA=False)
(-27*a*d*thetainv + (-3*b^2 - 9*a*c)*one - 9*a*b*theta)*x^2 + (-18*b*d*thetainv - 12*b*c*one - 18*a*c*theta)*x*y + (-9*c*d*thetainv + (-3*c^2 - 9*b*d)*one - 27*a*d*theta)*y^2
```

The FiniteDimensionalAlgebra class has a lot of useful functionality you
can get at using the algebra method of BinaryCubicForm. Like if you want
to quickly find out what the general minimal polynomials of &omega; and
&theta; are:
```
sage: one, omega, theta = BCF2.algebra().basis()
sage: omega.minimal_polynomial()
x^3 + b*x^2 + a*c*x + a^2*d
sage: theta.minimal_polynomial()
x^3 - c*x^2 + b*d*x - a*d^2
```

You can compute the number of subrings of the associated cubic ring of
given prime index (for forms over **Z**):
```
sage: BCF.number_of_index_p_subrings(2)
3
sage: BCF1.number_of_index_p_subrings(2)
1
```

And also the number of index *p* overrings:
```
BCF.number_of_p_overrings(2)
0
BCF1.number_of_p_overrings(3)
1
```

You can in fact get a list of primes at which the form is non-maximal
(or really its associated cubic ring):
```
sage: BCF.non_maximal_primes()
[]
sage: BCF1.non_maximal_primes()
[3]
```

## Quartic rings and pairs of ternary quadratic forms

In this section, I'll give some examples of using the Bhargavology code to study quartic rings. Again, I've tried to include some interesting examples. The mathematical details for the relation between quartic rings and pairs of ternary quadratic forms are mostly to be found in Manjul Bhargava's article [Higher composition laws III: The parametrization of quartic rings](https://doi.org/10.4007/annals.2004.159.1329).

### Some basic functionality

One way to create a pair of ternary quadratic forms is to pass it a list of coefficients. This list has the coefficients of one ternary quadratic form and then the other.
```
sage: R = PolynomialRing(QQ, names=['a%s%s'%(i, j) for i in range(3) for j in range(3)] + ['b%s%s'%(i, j) for i in range(3) for j in range(3)])
sage: R.inject_variables()
Defining a00, a01, a02, a10, a11, a12, a20, a21, a22, b00, b01, b02, b10, b11, b12, b20, b21, b22
sage: coeffs = [a00, a11, a22, a01, a12, a02, b00, b11, b22, b01, b12, b02]; coeffs
[a00, a11, a22, a01, a12, a02, b00, b11, b22, b01, b12, b02]
sage: PairOfTernaryQuadraticForms(coeffs, base_ring=R).print_pair()
[    a00 1/2*a01 1/2*a02|    b00 1/2*b01 1/2*b02]
[1/2*a01     a11 1/2*a12|1/2*b01     b11 1/2*b12]
[1/2*a02 1/2*a12     a22|1/2*b02 1/2*b12     b22]
```

The simplest way to create a pair of ternary quadratic forms is probably from a quartic polynomial:
```
sage: x = polygen(ZZ)
sage: f = x^4 - x^3 - 3*x^2 + 2*x + 4
sage: PTQF = _naive_pair_from_poly(f); PTQF
[-x1*x2 + x2^2 + 2*x1*x3 - 3*x2*x3 + 4*x3^2, -x1^2 + x2*x3]
```

This gives a pair of ternary quadratic forms corresponding to the quartic ring generated by f. We can see this by working with a general quartic. Specifically, we can create a quartic algebra with basis the normalized basis of the pair of ternary quadratic forms obtained from f, and we can see that the negative of the third basis element has characteristic polynomial f:
```
sage: R.<a0, a1, a2, a3> = PolynomialRing(QQ)
sage: Rx.<x> = PolynomialRing(R)
sage: f2 = Rx(list(R.gens()) + [1]); f2
x^4 + a3*x^3 + a2*x^2 + a1*x + a0
sage: PTQF2 = _naive_pair_from_poly(f2); PTQF2
[a3*x1*x2 + x2^2 + a1*x1*x3 + a2*x2*x3 + a0*x3^2, -x1^2 + x2*x3]
sage: A = PTQF2.algebra()
sage: (-A.basis()[2]).characteristic_polynomial()
x^4 + a3*x^3 + a2*x^2 + a1*x + a0
```

We can compare discriminants:
```
sage: f.discriminant()
a1^2*a2^2*a3^2 - 4*a0*a2^3*a3^2 - 4*a1^3*a3^3 + 18*a0*a1*a2*a3^3 - 27*a0^2*a3^4 - 4*a1^2*a2^3 + 16*a0*a2^4 + 18*a1^3*a2*a3 - 80*a0*a1*a2^2*a3 - 6*a0*a1^2*a3^2 + 144*a0^2*a2*a3^2 - 27*a1^4 + 144*a0*a1^2*a2 - 128*a0^2*a2^2 - 192*a0^2*a1*a3 + 256*a0^3
sage: PTQF2.disc()
a1^2*a2^2*a3^2 - 4*a0*a2^3*a3^2 - 4*a1^3*a3^3 + 18*a0*a1*a2*a3^3 - 27*a0^2*a3^4 - 4*a1^2*a2^3 + 16*a0*a2^4 + 18*a1^3*a2*a3 - 80*a0*a1*a2^2*a3 - 6*a0*a1^2*a3^2 + 144*a0^2*a2*a3^2 - 27*a1^4 + 144*a0*a1^2*a2 - 128*a0^2*a2^2 - 192*a0^2*a1*a3 + 256*a0^3
```

We can also display a pair of ternary quadratic forms as the pair of corresponding Gram matrices:
```
sage: PTQF2.print_pair()
[     0 1/2*a3 1/2*a1|    -1      0      0]
[1/2*a3      1 1/2*a2|     0      0    1/2]
[1/2*a1 1/2*a2     a0|     0    1/2      0]
```
(Using PTQF2.show() prints the pair out all pretty!)

We can ask for a multiplication table:
```
sage: PTQF2.multiplication_table()
                | 1             alpha1   alpha2   alpha3
+---------------+-------------+--------+--------+--------+
  alpha1^2      | -a1*a3 - a0   -a2      a1       -a3
  alpha1*alpha2 | -a1           0        0        -1
  alpha1*alpha3 | -a0*a3        0        a0       -a2
  alpha2^2      | -a2           -1       a3       0
  alpha2*alpha3 | -a0           0        0        0
  alpha3^2      | 0             a0       0        -a1
```

There are several invariants/covariants one has access to (for definitions see [HCLIII], or page 114 of Bhargava's PhD thesis for the quadratic covariant, which is a GL(3)-equivariant ternary quadratic form):
```
sage: PTQF2.lambda_ijkl([1,1,2,3])
a2
sage: PTQF2.cijk([1,1,3])
-a3
sage: PTQF2.cubic_invariant()
(a1*a2*a3 - a0*a3^2 - a1^2)*X^3 + (-a2^2 - a1*a3 + 4*a0)*X^2*Y + 2*a2*X*Y^2 - Y^3
sage: PTQF2.quadratic_covariant()
(4*a2^2 - 8*a1*a3 - 16*a0)*x1^2 + (4*a2*a3 - 24*a1)*x1*x2 + (3*a3^2 - 8*a2)*x2^2 + (4*a1*a2 - 24*a0*a3)*x1*x3 + (2*a1*a3 - 32*a0)*x2*x3 + (3*a1^2 - 8*a0*a2)*x3^2
```

We can act by GL(2), GL(3), or a pair (g3, g2) where g3 is in GL(3) and g2 is in GL(2):
```
sage: g2 = Matrix(2, [1, 5, 0, 2])
sage: g3 = Matrix(3, [2, 0, 1, 0, 3, 0, 1, 0, 0])
sage: PTQF2.act_by_GL2(g2).print_pair()
[    -6     a3     a1|    -2      0      0]
[    a3      2 a2 + 3|     0      0      1]
[    a1 a2 + 3   2*a0|     0      1      0]
sage: PTQF2.act_by_GL2(g3).print_pair()
[   0   a3   a1|  -3    0    0]
[  a3    2   a2|   0    0  3/2]
[  a1   a2 2*a0|   0  3/2    0]
sage: PTQF2.action([g3, g2]).print_pair()
[2*a0 + 4*a1 - 24  3*a2 + 6*a3 + 9          a1 - 12|              -8                3               -4]
[ 3*a2 + 6*a3 + 9               18             3*a3|               3                0                0]
[         a1 - 12             3*a3               -6|              -4                0               -2]
```

There are more functions available. Please see the source code if you're interested!
