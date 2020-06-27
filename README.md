# Bhargavology
Sage (python) scripts for working with cubic and quartic rings via binary cubic forms and pairs of ternary quadratic forms in the spirit of Delone–Faddeev, Davenport–Heilbronn, Gan–Gross–Savin, Bhargava, et al.

**Note:** this code is still in development. I'm releasing it because I think it can be useful to people out there despite this fact.

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

When a binary cubic forms corresponds to an order in a number field (i.e. when the form
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
