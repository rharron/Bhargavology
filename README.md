# Bhargavology
Sage (python) scripts for working with cubic and quartic rings via binary cubic forms and pairs of ternary quadratic forms in the spirit of Delone–Faddeev, Davenport–Heilbronn, Gan–Gross–Savin, Bhargava, et al.

**Note:** this code is still in development. I'm releasing it because I think it can be useful to people out there despite this fact.

## Initial usage comments

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

### Dedekind's example

One major advantage of using binary cubic forms is that it allows you to work with non-monogenic cubic rings. In this example, I'll use Dedekind's original example of a non-monogenic ring of integers to illustrate some functionality of this Bhargavology code. Dedekind's example is  **Q**(&alpha;) where &alpha; is a root of *f*(*x*)&nbsp;=&nbsp;*x*<sup>3</sup>&minus;*x*<sup>2</sup>&minus;2*x*&minus;8.

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
We could have also gotten this "by hand" as follows. Factoring the discriminant of BCF, we can see that the only square dividing it is 2<sup>2</sup>, so if *O* is not maximal, then it has index 2 in the maximal order. Because the last coefficient of BCF is divisible by 2<sup>2</sup> and the second-to-last by 2, finding a form corresponding to a ring containing O with index 2 is straightforward:

```
sage: BCF.disc().factor()
-1 * 2^2 * 503
sage: BCF.action(Matrix(2, [1,0,0,1/2]))
2*X^3 - X^2*Y - X*Y^2 - 2*Y^3
```

We can see the multiplication table *O* in the normalized basis 1, &omega;, &theta; (see section 2 of Bhargava–Shankar–Tsimerman, where they use the term "normal basis"):
```
sage: BCF.multiplication_table()
        | 1       omega               theta
+-------+-------+-------------------+-----------------------+
  1     | 1       omega               theta
  omega | omega   omega + theta + 2   8
  theta | theta   8                   8*omega - 2*theta - 8
```

Now that we have a binary cubic form corresponding to the maximal order, we can use it to understand some algebraic number theory. The rest of this example doesn't explicitly use the Bhargavology code, but I wanted to give an idea of what you can do with binary cubic forms. Dedekind's example is one where 2 always divides the index of any monogenic order, but using binary cubic forms, you can still determine the prime factorization of 2 in a manner analogous to Dedekind's theorem for polynomials (see Section 3 of del Corso–Dvornicich–Simon's "Decomposition of primes in non-maximal orders" for the mathematical details). Specifically, factoring BCFmax mod 2 gives polynomials that you should use to generate the prime ideals dividing 2*O*<sub>max</sub>:

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
