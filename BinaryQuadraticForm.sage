"""
@author: Robert Harron

Sage (python) script for working with binary quadratic forms as
they relate to binary cubic forms. Also contains the function
``Siegel_reduced_forms`` for enumerating indefinite binary
quadratic forms over ZZ that are reduced in the sense of Siegel.

"""

class BinaryQuadraticForm(SageObject):
    r"""
    A class representing a binary quadratic form over a ring R.

    EXAMPLES::

        Creating a binary quadratic form over ZZ is as simple as passing a list of coefficients.

        sage: BQFZZ = BinaryQuadraticForm([1, -1, 1]); BCFZZ
        x^2 - x*y + y^2

        One can create binary quadratic forms over other rings as well, and use different variable names.

        sage: K.<sq2> = QuadraticField(2)
        sage: R = K.maximal_order()
        sage: BinaryCubicForm([1, sq2, 1 + sq2], base_ring=R, var_names='X, Y')
        ^2 + sq2*X*Y + (sq2 + 1)*Y^2

        One can also work with a generic binary quadratic form, and also act by invertible 2x2 matrices.

        sage: R.<a, b, c> = PolynomialRing(QQ)
        sage: BQF = BinaryQuadraticForm([a, b, c], base_ring=R); BQF
        a*x^2 + b*x*y + c*y^2
        sage: BQF.disc()
        b^2 - 4*a*c
        sage: S.<A,B,C,D> = PolynomialRing(R)
        sage: BQF.action(Matrix(2, [A, B, C, D]))
        (a*A^2 + b*A*B + c*B^2)*x^2 + (2*a*A*C + b*B*C + b*A*D + 2*c*B*D)*x*y + (a*C^2 + b*C*D + c*D^2)*y^2

    ..TODO::

        Add more examples!
    """
    def __init__(self, coeffs, base_ring=ZZ, var_names='x, y', adjuster=None):
        r"""
        Initializes this instance of the BinaryQuadraticForm class.

        INPUT:

        - ``coeffs`` -- a listable with 3 entries giving the coefficients. The ordering is such that
        [a, b, c] corresponds to the binary quadratic form aX^2 + bXY + cY^2.
        - ``base_ring`` -- (Default: ZZ) the ring of coefficients for this binary quadratic form.
        - ``var_names`` -- (Default: 'x, y') the variable names to use for this binary quadratic form
        - ``adjuster`` -- (Default: None) If None, nothing special happens. Otherwise this should be a
        function that takes in a 2x2 invertible matrix and applies an automorphism to it. This is
        applied to a matrix before it acts on this form.
        """
        self._a, self._b, self._c = [base_ring(c) for c in coeffs]
        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(base_ring, var_names)
        self._vars = self._polynomial_ring.gens()
        X, Y = self._vars
        self._polynomial = self.a * X^2 + self.b * X * Y + self.c * Y^2
        self._names = var_names
        self._adjuster = adjuster
        self._UHP = HyperbolicPlane().UHP() # Only some times used. TODO: move this out of this class
        return

    def __call__(self, x, y):
        r"""
        Evaluates the binary quadratic form at (x, y).

        INPUT:

        - ``x`` -- something Sage's coercion system can plug into the polynomial of ``self``
        - ``y`` -- something Sage's coercion system can plug into the polynomial of ``self``

        OUTPUT:

        This binary quadratic form evaluated at the point (x, y)
        """
        return self._polynomial(x, y)

    def base_ring(self):
        r"""
        Return the base ring of this binary quadratic form, i.e. the ring its coefficients are elements of.
        """
        return self._base_ring

    def change_ring(self, new_base):
        r"""
        Change the base ring of this binary quadratic form.

        INPUT:

        - ``new_base`` -- a ring to which Sage can convert the current coefficients of ``self``.
        """
        return BinaryQuadraticForm([self._a, self._b, self._c], base_ring=new_base)

    def disc(self):
        r"""
        Return the discriminant of this binary quadratic form. The output will be an element of the base ring of ``self``.
        """
        return self.b^2 - 4 * self.a * self.c

    def coefficients(self):
        r"""
        Return the coefficients of ``self`` as a tuple. The order is the same as that used in ``__init__``.
        """
        return self.a, self.b, self.c

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

    def polynomial(self):
        r"""
        Return the underlying polynomial of this binary quadratic form.
        """
        return self._polynomial

    def action(self, gamma):
        r"""
        Return the binary quadratic form obtained by acting on ``self`` via a linear change
        of variables given by applying ``self._adjuster`` to the invertible 2x2 matrix
        ``gamma``.

        Specifically, if F(x,y) denotes this binary quadratic form, then this function
        returns F((x,y) * g' ), where g' = ``self._adjuster(gamma)`` and (x,y) * g' denotes
        row vector-matrix multiplication.

        INPUT:

        - ``gamma`` -- an invertiable 2x2 matrix.

        ..TODO::

        - Use correct variable names.
        - Use the coercion framework better.
        - Try to return a better base ring.
        """
        X, Y = self._vars
        S = gamma[0,0].parent()
        R = X.parent().change_ring(S)
        if not self._adjuster is None:
            gamma = self._adjuster(gamma)
        Xp, Yp = vector(R, [X, Y]) * gamma
        Fp = self._polynomial(Xp, Yp)
        X, Y = Fp.variables()
        return BinaryQuadraticForm([Fp[X^(2-i) * Y^i] for i in range(3)], base_ring=S)

    def gram_matrix(self):
        r"""
        Return the gram_matrix associated to ``self``.

        Try to return it over ``self.base_ring()``; otherwise use the fraction
        field.
        """
        try:
            return Matrix(self._base_ring, 2, [self.a, self.b/2, self.b/2, self.c])
        except(TypeError):
            return Matrix(self._base_ring.fraction_field(), 2, [self.a, self.b/2, self.b/2, self.c])

    def geodesic(self, ambient=QQbar):
        r"""
        Return the geodesic attached to this binary quadratic form.
        The form must be an indefinite form over a ring that can
        convert into ``ambient``.

        INPUT:

        - ``ambient`` -- (Default: QQbar) a ring that contains the roots
        in PP^1 of this binary quadratic form.
        """
        Rx = PolynomialRing(self._base_ring, 'x')
        x, y = self._vars
        P = self._polynomial
        rs = []
        if y.divides(P):
            rs.append(Infinity)
        P = Rx(P(x, 1))
        Prs = P.roots(ambient)
        rs += [r for r, _ in Prs]
        return self._UHP.get_geodesic(*rs)

    def circle(self, color='blue', fontsize=10, rotation=0.0):
        r"""
        Return a graphics object that is the geodesic in the upper-half
        plane corresponding to ``self``. The form must be an indefinite form.
        """
        centre = -self.b / (2 * self.a)
        r = self.disc().sqrt() / (2 * self.a.abs())
        return arc((centre, 0), r, sector=(0,pi), color=color)#+text('(%s,%s,%s)'%(self._a, self._b, self._c),(centre,r + sign(self._a)*0.01*fontsize), fontsize=fontsize, color=color, rotation=rotation)

    def radius(self):
        r"""
        Return the radius of the geodesic in the upper-half plane corresponding
        to ``self``. The form must be an indefinite form.
        """
        return self.disc().sqrt() / (2 * self.a.abs())

    def centre(self):
        r"""
        Return the centre of the geodesic in the upper-half plane corresponding
        to ``self``. The form must be an indefinite form.
        """
        return -self.b / (2 * self.a)

    def zvalue(self):
        r"""
        Return the point in the upper-half plane corresponding to this
        binary quadratic form. The form must be positive definite. The output
        lies in QQbar.
        """
        i = QQbar(-1).sqrt()
        a, b, c = [QQbar(cc) for cc in self.coefficients()]
        return (b + i * AA(4 * a * c - b^2).sqrt()) / (2 * a)

def Siegel_reduced_forms(D):
    r"""
    Return a list of triples (a,b,c) giving the coefficients of all Siegel-reduced
    indefinite binary quadratic forms of discriminant ``D``.

    An indefinite binary quadratic form over ZZ is Siegel-reduced if its
    corresponding geodesic intersects the usual Gauss fundamental domain. See
    Siegel's Lecture on Quadratic Forms for details.

    INPUT:

    - ``D`` -- a positive discriminant

    EXAMPLES::

        sage: Siegel_reduced_forms(5)
        [(1, 1, -1), (-1, 1, 1), (1, -1, -1), (-1, -1, 1)]
    """
    reds = []
    D = ZZ(D)
    abound = (D/3).sqrt().n().floor()
    for b in srange((4*D/3).sqrt().n().floor() + 1):
        for a in srange(1, abound + 1):
            c = (b^2 - D) / (4 * a)
            if c not in ZZ:
                continue
            c = ZZ(c)
            if 3*a*c > D/4:
                print a,b,c
                continue
            reds.append((a,b,c))
            reds.append((-a,b,-c))
            if b != 0:
                reds.append((a,-b,c))
                reds.append((-a,-b,-c))
    return [(a,b,c) for a,b,c in reds if ((1/2 + b/(2*a))^2+3/4 <= D/(4*a^2)) or ((-1/2 + b/(2*a))^2+3/4 <= D/(4*a^2))]
