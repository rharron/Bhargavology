"""
@author: Robert Harron

Sage (python) script for working with binary cubic forms and
cubic rings.

"""

load('BinaryQuadraticForm.sage')

from sage.categories.magmatic_algebras import MagmaticAlgebras

class BinaryCubicForm(SageObject):
    r"""
    A class representing a binary cubic form over a ring R.

    EXAMPLES::

        Creating a binary cubic form over ZZ is as simple as passing a list of coefficients.

        sage: BCFZZ = BinaryCubicForm([1, -1, -2, -8]); BCFZZ
        X^3 - X^2*Y - 2*X*Y^2 - 8*Y^3

        One can create binary cubic forms over other rings as well, and use different variable names.

        sage: K.<I> = QuadraticField(-1)
        sage: R = K.maximal_order()
        sage: BinaryCubicForm([1, -I, 2+I, -8], base_ring=R, var_names='x_1, x_2')
        x_1^3 + (-I)*x_1^2*x_2 + (I + 2)*x_1*x_2^2 + (-8)*x_2^3

        One can also work with a generic binary cubic form.

        sage: R.<a, b, c, d> = PolynomialRing(QQ)
        sage: BCF = BinaryCubicForm([a, b, c, d], base_ring=R); BCF
        a*X^3 + b*X^2*Y + c*X*Y^2 + d*Y^3
        sage: BCF.disc()
        b^2*c^2 - 4*a*c^3 - 4*b^3*d + 18*a*b*c*d - 27*a^2*d^2
        sage: BCF.hessian()
        (b^2 - 3*a*c)*X^2 + (b*c - 9*a*d)*X*Y + (c^2 - 3*b*d)*Y^2
        sage: BCF.multiplication_table()
                | 1       omega                        theta
        +-------+-------+----------------------------+----------------------------+
          1     | 1       omega                        theta
          omega | omega   (-b)*omega + a*theta - a*c   -a*d
          theta | theta   -a*d                         (-d)*omega + c*theta - b*d

        The ``action`` method gives the twisted action by GL(2).
        sage: R2.<N> = PolynomialRing(R)
        sage: BCF.action(Matrix(2, [1,N,0,1]))
        (d*N^3 + c*N^2 + b*N + a)*X^3 + (3*d*N^2 + 2*c*N + b)*X^2*Y + (3*d*N + c)*X*Y^2 + d*Y^3
        sage: BCFZZ.action(Matrix(2, [1,0,0,1/2]))
        2*X^3 - X^2*Y - X*Y^2 - 2*Y^3

    ..TODO::

        Add more examples!
    """

    def __init__(self, coeffs, base_ring=ZZ, var_names='X, Y'):
        r"""
        Initializes this instance of the BinaryCubicForm class.

        INPUT:

        - ``coeffs`` -- a listable with 4 entries giving the coefficients. The ordering is such that
        [a, b, c, d] corresponds to the binary cubic form aX^3 + bX^2Y + cXY^2 + dY^3.
        - ``base_ring`` -- (Default: ZZ) the ring of coefficients for this binary cubic form. A binary
        cubic form over a ring R will correspond to a cubic algebra over R.
        - ``var_names`` -- the variable names to use for this binary cubic form
        """
        self._a, self._b, self._c, self._d = [base_ring(c) for c in coeffs]
        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(base_ring, var_names)
        self._vars = self._polynomial_ring.gens()
        self._var_names = var_names
        X, Y = self._vars
        self._polynomial = self._a * X^3 + self._b * X^2 * Y + self._c * X * Y^2 + self._d * Y^3
        return

    def __call__(self, x, y):
        r"""
        Evaluates the binary cubic form at (x, y).

        INPUT:

        - ``x`` -- something Sage's coercion system can plug into the polynomial of ``self``
        - ``y`` -- something Sage's coercion system can plug into the polynomial of ``self``

        OUTPUT:

        This binary cubic form evaluated at the point (x, y)
        """
        return self._polynomial(x, y)

    def base_ring(self):
        r"""
        Return the base ring of this binary cubic form, i.e. the ring its coefficients are elements of.
        """
        return self._base_ring

    def change_ring(self, new_base):
        r"""
        Change the base ring of this binary cubic form.

        INPUT:

        - ``new_base`` -- a ring to which Sage can convert the current coefficients of ``self``.
        """
        return BinaryCubicForm([self._a, self._b, self._c, self._d], base_ring=new_base)

    def disc(self):
        r"""
        Return the discriminant of this binary cubic form. The output will be an element of the base ring of ``self``.
        """
        return self._b^2*self._c^2-4*self._a*self._c^3-4*self._b^3*self._d-27*self._a^2*self._d^2+18*self._a*self._b*self._c*self._d

    def coefficients(self):
        r"""
        Return the coefficients of ``self`` as a tuple. The order is the same as that used in ``__init__``.
        """
        return self._a, self._b, self._c, self._d

    def polynomial(self):
        r"""
        Return the underlying polynomial of this binary cubic form.
        """
        return self._polynomial

    def hessian(self, var_names=None):
        r"""
        Return the Hessian of this binary cubic form as an object of class BinaryQuadraticForm.

        The Hessian a binary cubic form is the determinant of the matrix of its second order partial derivatives. This is a classical
        quadratic covariant of the binary cubic form. In more modern terms, if the space of binary cubic forms is given the twisted
        action of GL(2) as in Bhargava–Shankar–Tsimerman and the space of binary quadratic forms is given the usual action, then the
        Hessian is GL(2)-equivariant.

        INPUT:

        - ``var_names`` -- (Default: None) the variable names to use for the Hessian. If None, the same variables as ``self`` will be used.
        """
        if var_names is None:
            var_names = ','.join(self._polynomial_ring.variable_names())
        return BinaryQuadraticForm([self._b^2-3*self._a*self._c, self._b*self._c-9*self._a*self._d, self._c^2-3*self._b*self._d], base_ring=self._base_ring, var_names=var_names)

    def shape(self, check=True, names=["thetainv", "one", "theta"], use_AA=True):
        r"""
        Return (a binary quadratic form representing) the shape of this binary cubic form.

        This is primarily meant to be for irreducible binary cubic forms over ZZ or QQ,
        in which case the shape is the binary quadratic form representing the Minkowski
        inner product with respect to the normalized basis of ``self``. There is a formula
        for this binary quadratic form in Chapter 9 of Terr's thesis; it is given in terms
        of the coefficients of ``self`` and a root of the associated polynomial. In other
        cases, this formula is used as a definition for what is returned.

        INPUT:

        - ``check`` -- (Default: True) whether to check if the discriminant of ``self`` is
        positive. Designed to be used only when ``self`` is over ZZ or QQ.
        - ``names`` -- (Default: ["thetainv", "one", "theta"])
        only used when use_AA is False and either ``check`` is ``False`` or ``self.disc()`` is
        non-positive. In this case, this function returns a binary quadratic form whose
        base ring is a ``FiniteDimensionalAlgebra`` over the field of fractions of
        ``self.base_ring()``. The parameter ``names`` is then used as the display names
        of the basis of this finite-dimensional algebra. In this basis, the middle variable
        is the identity, the product of the other two is the identity, and the last
        variable is a root of F(x, 1), where F is this binary cubic form.
        -use_AA=True -- (Default: True) ignored if ``check`` is ``True`` and ``self.disc()``
        is positive. Meant to be used when ``self.base_ring()`` is ZZ or QQ, or any ring
        in a number field. If True, a root of F(x, 1) is found in AA and the shape is
        returned as a binary quadratic form over AA. If false, the finite-dimensional
        algebra method described in the ``names`` parameter is used.

        OUTPUT:

        If ``check`` is ``True`` and ``self.disc()`` is positive, then this returns
        ``self.hessian()``.

        Otherwise, if ``use_AA`` is ``True``, then a root of F(x, 1) is found in AA and
        the shape of ``self`` is returned as a ``BinaryQuadraticForm`` over AA.

        Otherwise, ``self.base_ring()`` must be an integral domain R. A finite-dimensional
        algebra A over the fraction field of R containing a root theta of F(x, 1) is then
        constructed and the shape is returned as a BinaryQuadraticForm over A.

        .. WARNING::

            This function has not really been tested beyond its use cases.
        """
        #Fix this maybe for reducible and/or degenerate forms
        if check and self.disc() > 0:
            return self.hessian()
        a, b, c, d = self.coefficients()
        if use_AA:
            x = polygen(QQ)
            theta = self._polynomial(x, 1).roots(AA)[0][0]
            thetainv = ~theta
            one = 1
            return BinaryQuadraticForm([-3*(9*a*d*(thetainv)+(b^2+3*a*c)*one+3*a*b*theta), -6*(3*b*d*(thetainv)+2*b*c*one+3*a*c*theta), -3*(3*c*d*(thetainv)+(c^2+3*b*d)*one+9*a*d*theta)], base_ring=AA)
        K = self._base_ring.fraction_field()
        matrices = [Matrix(K, 3, [-c/d, -b/d, -a/d, 1,0,0, 0,1,0]), identity_matrix(3), Matrix(K, 3, [0,1,0, 0,0,1, -d/a, -c/a, -b/a])]
        A = FiniteDimensionalAlgebra(K, matrices, assume_associative=True, names=names, category=MagmaticAlgebras(K).FiniteDimensional().WithBasis().Associative().Commutative().Unital().Inverse())
        thetainv, one, theta = A.basis()
        return BinaryQuadraticForm([-3*(9*a*d*(thetainv)+(b^2+3*a*c)*one+3*a*b*theta), -6*(3*b*d*(thetainv)+2*b*c*one+3*a*c*theta), -3*(3*c*d*(thetainv)+(c^2+3*b*d)*one+9*a*d*theta)], base_ring=A)

    def trace_zero_form(self):
        r"""
        Return a BinaryQuadraticForm representing the trace zero form of the
        cubic ring associated to ``self`` with respect to its normalized basis.
        This is 6 times the Hessian.
        """
        return BinaryQuadraticForm([6*(self._b^2-3*self._a*self._c), 6*(self._b*self._c-9*self._a*self._d), 6*(self._c^2-3*self._b*self._d)], self._base_ring)

    def nice_basis(self, names=['omega', 'theta']):
        r"""
        Return a normalized basis of ``self``, where ``self`` should have
        base ring ZZ or QQ and be irreducible.

        INPUT:

        - ``names`` -- (Default: ['omega', 'theta']) only the first component
        of ``names`` is used. It is passed to NumberField as the name for the
        generator.

        OUTPUT:

        Writing F for this binary cubic form, this function creates a NumberField K
        whose generator is named ``names[0]``. The function returns the normalized
        basis (1, omega, theta) of the order in K corresponding to F, where omega is
        the generator of K.
        """
        x = polygen(self._base_ring)
        a = self._a
        f = a^2 * self._polynomial(x / a, 1)
        K = NumberField(f, names[0])
        omega = K.gen()
        eta = omega / a
        theta = -self._d / eta
        return (K.one(), omega, theta)

    def algebra(self, names=['one', 'omega', 'theta']):
        r"""
        Returns a FiniteDimensionalAlgebra over the fraction field of ``self.base_ring()``
        (so the latter must be an integral domain). This algebra has basis the normalized
        basis of the cubic ring corresponding to ``self``.

        INPUT:

        - ``names`` -- (Default: ['one', 'omega', 'theta']) the names to assign to the
        basis elements of the outputted algebra.
        """
        a, b, c, d = self.coefficients()
        m = -a*c
        n = -a*d
        ell = -b*d
        matrices = [identity_matrix(3), matrix(3, [0,1,0, m,-b,a, n,0,0]), matrix(3, [0,0,1, n,0,0, ell,-d,c])]
        #for m in matrices:
        #    show(m)
        return FiniteDimensionalAlgebra(self._base_ring.fraction_field(), matrices, assume_associative=True, names=names)

    def order(self, with_basis=False):
        r"""
        Returns the cubic order corresponding to ``self``. In particulat, ``self``
        should be over ZZ or QQ and be irreducible.

        INPUT:

        - ``with_basis`` -- (Default: False) if true, return a tuple where the first
        element is the order usually returned and the second element is itself a tuple
        (omega, theta), where 1, omega, theta form a normalized basis of the order.

        OUTPUT:

        Letting F denote this binary cubic form, this function creates a number field K
        generated by an element eta that is a root of F(x, 1). The order it returns is
        the order in K generated by omega and theta, where omega and theta give the
        normalized basis of the order corresponding to F. Specifically, omega = a*eta
        and theta = -d / eta.
        """
        x = polygen(QQ)
        f = self.polynomial()(x, 1)
        K.<eta> = NumberField(f)
        omega = self._a * eta
        theta = -self._d / eta
        O = K.order([omega, theta])
        if with_basis:
            return O, (omega, theta)
        return O

    def _repr_(self):
        r"""
        Returns the string representation of polynomial underlying ``self``.
        """
        return str(self._polynomial)

    def _latex_(self):
        r"""
        Returns the latex string representation of polynomial underlying ``self``.
        """
        return self._polynomial._latex_()

    def multiplication_table(self, names=['omega', 'theta']):
        r"""
        Return a Sage ``table`` displaying the multiplication table of the normalized
        basis of ``self``.

        INPUT:

        - ``names`` -- (Default: ['omega', 'theta']) the names to display for the
        non-identity elements of the normalized basis.
        """
        a, b, c, d = self.coefficients()
        n = -a * d
        m = -a * c
        ell = -b * d
        R = PolynomialRing(self.base_ring(), names)
        omega, theta = R.gens()
        return table([['', '1'] + names, ['1', '1'] + names, [names[0], names[0], str(m - b*omega + a*theta), str(n)], [names[1], names[1], str(n), str(ell - d*omega + c*theta)]], header_row=True, header_column=True)

    def twisted_action(self, gamma):
        r"""
        Return the binary cubic form obtained by acting on ``self`` via the so-called
        twisted left action of the invertible 2x2 matrix ``gamma``.

        Specifically, if F(x,y) denotes this binary cubic form and g denotes ``gamma``, then
        this function returns F((x,y) *g ) / det(g), where (x,y) * g denotes
        row vector-matrix multiplication. This function attempts to return a binary cubic form
        whose base ring is the same as ``self``, otherwise it tries to use the base ring of
        ``gamma``.

        INPUT:

        - ``gamma`` -- an invertiable 2x2 matrix.

        ..TODO::

        - Use correct variable names.
        - Use the coercion framework better.
        """
        X, Y = self._vars
        S = gamma.base_ring()
        R = X.parent().change_ring(S)
        Xp, Yp = vector(R, [X, Y]) * gamma
        Fp = self._polynomial(Xp, Yp) / gamma.det()
        X, Y = Fp.variables()
        try:    #First try keeping self's base ring
            return BinaryCubicForm([Fp[X^(3-i) * Y^i] for i in range(4)], base_ring=self._base_ring)
        except TypeError:
            pass
        try:    #Then try using gamma's base ring
            return BinaryCubicForm([Fp[X^(3-i) * Y^i] for i in range(4)], base_ring=S)
        except TypeError:
            return BinaryCubicForm([Fp[X^(3-i) * Y^i] for i in range(4)], base_ring=S.fraction_field())

    def action(self, gamma):
        r"""
        Returns ``self.twisted_action(gamma)``.

        INPUT:

        - ``gamma`` -- an invertiable 2x2 matrix.

        ..TODO::

        -Implement an "adjuster" as in BinaryQuadraticForm
        """
        return self.twisted_action(gamma)

    def untwisted_action(self, gamma):
        r"""
        Return the binary cubic form obtained by acting on ``self`` via the usual
        left action of the invertible 2x2 matrix ``gamma``.

        Specifically, if F(x,y) denotes this binary cubic form and g denotes ``gamma``, then
        this function returns F((x,y) *g ), where (x,y) * g denotes
        row vector-matrix multiplication.

        INPUT:

        - ``gamma`` -- an invertiable 2x2 matrix.

        ..TODO::

        - Try to return a better base ring as in ``twisted_action``.
        - Use correct variable names.
        - Use the coercion framework better.
        """
        X, Y = self._vars
        S = gamma[0,0].parent()
        R = X.parent().change_ring(S)
        Xp, Yp = vector(R, [X, Y]) * gamma
        Fp = self._polynomial(Xp, Yp)
        X, Y = Fp.variables()
        return BinaryCubicForm([Fp[X^(3-i) * Y^i] for i in range(4)], base_ring=S)

    def number_of_index_p_subrings(self, p):
        r"""
        Return the number of index ``p`` subrings in the cubic ring corresponding to
        this binary cubic form, where ``p`` is a prime number.

        This function requires ``self.base_ring()`` to be ZZ. The value returned is
        taken from Proposition 15 of BST.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("number_of_index_p_subrings currently implemented only for binary cubic forms over ZZ.")
        if all([c % p == 0 for c in self.coefficients()]):    #if self is 0 mod p
            return p + 1
        return sum([1 for F, _ in self.polynomial().change_ring(GF(p)).factor() if F.degree() == 1])

    def number_of_p_overrings(self, p):
        r"""
        Return the number of overrings of the cubic ring R corresponding to
        this binary cubic form that contain R with index ``p``, where ``p``
        is a prime number.

        This function requires ``self.base_ring()`` to be ZZ. The value returned is
        taken from Proposition 16 of BST.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("number_of_p_overrings currently implemented only for binary cubic forms over ZZ.")
        #First deal with case where self is 0 mod p
        if all([c % p == 0 for c in self.coefficients()]):
            count = 0
            good_root = True
            #First check the infinite root
            for i in range(p):
                if self(1, i*p) % p^2 != 0:
                    good_root = False
                    break
            if good_root:
                count += 1
            #Now check the non-infinite roots
            for j in range(p):
                good_root = True
                for i in range(p):
                    if self(j+i*p, 1) % p^2 != 0:
                        good_root = False
                        break
                if good_root:
                    count += 1
            return count

        #Now deal with self not 0 mod p
        X = polygen(GF(p))
        count = 0
        for F, e in self.polynomial().change_ring(GF(p)).factor():
            if e == 1:
                continue
            Y = self._vars[1].change_ring(GF(p))
            if F == Y:    #Deal with the root at infinity separately, if it exists
                if all(map(lambda x : (self(1, x) % p^2) == 0, [i * p for i in range(p)])):
                    count += 1
                continue
            f = F(X, 1)
            for r, _ in f.roots():    #Deal with the non-infinite roots
                if all(map(lambda x : (self(x, 1) % p^2) == 0, [r.lift_centered() + i * p for i in range(p)])):
                    count += 1
        return count

    def is_maximal(self, p=None):
        r"""
        Returns whether this binary cubic form corresponds to a maximal cubic ring, or
        if ``p`` is a prime number, whether it is p-maximal.

        This function requires ``self.base_ring()`` to be ZZ.

        INPUT:

        - ``p`` -- (Default: ``None``) ``None`` or a prime number.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("is_maximal currently implemented only for binary cubic forms over ZZ.")
        if p is not None:
            return self.number_of_p_overrings(p) == 0
        for p in self.disc().prime_divisors():
            if self.number_of_p_overrings(p) > 0:
                return False
        return True

    def non_maximal_primes(self):
        r"""
        Returns the prime numbers p for which the cubic ring corresponding to
        this binary cubic is not p-maximal.

        This function requires ``self.base_ring()`` to be ZZ.

        OUTPUT:

        A list of primes numbers.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("non_maximal_primes currently implemented only for binary cubic forms over ZZ.")
        return [p for p in self.disc().prime_divisors() if self.number_of_p_overrings(p) > 0]

    def find_p_overring_form(self, p):
        r"""
        Returns an overring of the cubic ring corresponding to ``self`` such that the latter
        ring has index ``p` in the former ring, where ``p`` is a prime number.

        This function requires ``self.base_ring()`` to be ZZ.

        OUTPUT:

        A BinaryCubicForm.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("find_p_overring_form currently implemented only for binary cubic forms over ZZ.")
        a, b, c, d = self.coefficients()
        if a % p^2 == 0 and b % p == 0:
            return self.superform(p)
        if d % p^2 == 0 and c % p == 0:
            return self.superform(p, reverse=True)
        #Deal with cases corresponding to non-zero, non-infinite roots
        #First deal with the case self is 0 mod p
        found_good_root = False
        if all([c % p == 0 for c in self.coefficients()]):
            for r in range(1, p):    #Only need to check non-zero, non-infinite roots
                good_root = True
                for i in range(p):
                    if self(r+i*p, 1) % p^2 != 0:
                        good_root = False
                        break
                if good_root:
                    found_good_root = True
                    break
            r = GF(p)(r)
        #Now deal with self not 0 mod p
        else:
            X = polygen(GF(p))
            for F, e in self.polynomial().change_ring(GF(p)).factor():
                if e == 1:
                    continue
                f = F(X, 1)
                for r, _ in f.roots():
                    if all(map(lambda x : (self(x, 1) % p^2) == 0, [r.lift_centered() + i * p for i in range(p)])):
                        found_good_root = True
                        break
                if found_good_root:
                    break
        if not found_good_root:
            raise ValueError("find_p_overring_form was given a p-maximal form (p = %s)"%p)
        rinv = (~r).lift_centered()
        BCFnew = self.action(Matrix(2, [1, rinv, 0, 1]))
        return BCFnew.superform(p)

    def superform(self, n, reverse=False):
        r"""
        Return a BinaryCubicForm whose corresponding normalized basis generates
        a lattice contaning the one generated by self's with index n in a very
        specific way, namely via the action of Matrix(2, [1,0,0,1/n]) (or
        if ``reverse`` is ``False``, then via the action of Matrix(2, [1/n,0,0,1])).
        """
        if reverse:
            return BinaryCubicForm([self._a * n, self._b, self._base_ring(self._c / n), self._base_ring(self._d / n^2)], base_ring=self._base_ring, var_names=self._var_names)
        return BinaryCubicForm([self._base_ring(self._a / n^2), self._base_ring(self._b / n), self._c, self._d * n], base_ring=self._base_ring, var_names=self._var_names)

    def find_maximal_form(self):
        r"""
        Return a binary cubic form GL(2,QQ)-equivalent to ``self`` that corresponds
        to a maximal cubic ring containing the ring corresponding to ``self``.

        This function requires ``self.base_ring()`` to be ZZ.

        ALGORITHM:

        This function iteratively uses ``find_p_overring_form(p)`` with the primes
        obtained from ``non_maximal_primes()`` until no more non maximal primes
        exist.

        OUTPUT:

        A BinaryCubicForm.
        """
        if not (self.base_ring() is ZZ):
            raise NotImplementedError("find_maximal_form currently implemented only for binary cubic forms over ZZ.")
        BCF = self
        bad_primes = BCF.non_maximal_primes()
        while len(bad_primes) > 0:
            for p in bad_primes:
                BCF = BCF.find_p_overring_form(p)
            bad_primes = BCF.non_maximal_primes()
        return BCF
