"""
@author: Robert Harron

Sage (python) script for working with ternary quadratic forms
as they relate to pairs of ternary quadratic forms.

Note: this is only meant to be used with PairOfTernaryQuadraticForms.sage.
It in fact references some functions defined in that file. So to use the
functionality in here simply do:
load('PairOfTernaryQuadraticForms.sage')
since that file loads this one.

"""

class TernaryQuadraticForm(SageObject):
    r"""
    A class for computing with ternary quadratic forms as they
    relate to the code from PairOfTernaryQuadraticForms.

    ..NOTE::This is currently a very bare bones implementation.
    """
    def __init__(self, coeffs, base_ring=None, var_names='x1, x2, x3'):
        r"""
        Initializes this instance of the TernaryQuadraticForm class.

        INPUT:

        - ``coeffs`` -- a listable with 6 entries giving the coefficients. The ordering is such that
        [a, b, c, d, e, f] corresponds to the ternary quadratic form
        a*x1^2 + d*x1*x2 + b*x2^2 + f*x1*x3 + e*x2*x3 + c*x3^2
        - ``base_ring`` -- (Default: None) the ring of coefficients for this ternary quadratic form.
        If ``None``, base_ring will be computed by using that of vector(coeffs).
        - ``var_names`` -- the variable names to use for this ternary quadratic form.

        EXAMPLES::

            sage: R.<a,b,c,d,e,f> = PolynomialRing(QQ)
            sage: Q = TernaryQuadraticForm([a,b,c,d,e,f], base_ring=R); Q
            a*x1^2 + d*x1*x2 + b*x2^2 + f*x1*x3 + e*x2*x3 + c*x3^2
        """
        if base_ring is None:
            base_ring = vector(coeffs).base_ring()
        self.coeffs = [base_ring(a) for a in coeffs]
        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(base_ring, var_names)
        self._vars = self._polynomial_ring.gens()
        self._polynomial = sum([self.coeffs[(-i^2 / 2 + 7 * i / 2) + j] * self._vars[j] * self._vars[j+i] for i in range(3) for j in range(3 - i)])
        return

    def _repr_(self):
        r"""
        Return the string representation of the underlying degree 2
        homogeneous polynomial of ``self``.

        EXAMPLES::

            sage: Q = TernaryQuadraticForm([1,2,3,4,5,6]); Q
            x1^2 + 4*x1*x2 + 2*x2^2 + 6*x1*x3 + 5*x2*x3 + 3*x3^2
        """
        return str(self._polynomial)

    def coefficients(self):
        r"""
        Return the coefficients of ``self`` as a list.

        EXAMPLES::

            sage: Q = TernaryQuadraticForm([1,2,3,4,5,6]); Q
            x1^2 + 4*x1*x2 + 2*x2^2 + 6*x1*x3 + 5*x2*x3 + 3*x3^2
            sage: Q.coefficients()
            [1, 2, 3, 4, 5, 6]
        """
        return self.coeffs

    def matrix(self):
        r"""
        Return the Gram matrix of ``self``.
        """
        try:
            A = _TQF_coeffs_to_matrix(self.coeffs, base_ring=self._base_ring)
        except TypeError:
            A = _TQF_coeffs_to_matrix(self.coeffs)
        return A

    def disc(self):
        r"""
        Return the discriminant of ``self``, i.e. 4 times the determinant
        of its Gram matrix.
        """
        return 4 * self.matrix().det()

    def action(self, gamma):
        r"""
        Return the TernaryQuadraticForm obtained by acting on ``self``
        by the 3x3 invertiable matrix ``gamma``.

        If A denotes the Gram matrix corresponding to ``self``, the
        Gram matric of the output is gamma * A * gamma^T; so this is
        the left action by GL(3).
        """
        A = self.matrix()
        gAgt = gamma * A * gamma.transpose()
        try:
            a2 = _matrix_to_TQF_coeffs(gAgt, base_ring=self._base_ring)
        except TypeError:
            a2 = _matrix_to_TQF_coeffs(gAgt)
        return TernaryQuadraticForm(a2, var_names=','.join(self._polynomial_ring.variable_names()))

    def primitive_form(self):
        r"""
        Return the primitive ternary quadratic form associated to ``self``, i.e. divide the
        coefficients by their gcd.
        """
        g = gcd(self.coeffs)
        return TernaryQuadraticForm([a/g for a in self.coeffs], base_ring=self._base_ring, var_names=','.join(self._polynomial_ring.variable_names()))
