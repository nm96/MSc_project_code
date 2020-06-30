"""Module containing functions which encode various nonlinear forcing models
which can be incorporated into the Jeffcott equations."""

def simple(X,dummy):
    """Model function for the trivial no-forcing model. Dummy parameter
    included in order to fit the convention for model-passing.
    Returns the tuple (fn,ft) = (normal force, tangential force).
    """
    return (0,0)

def VdH(X,h,k_c):
    """Model function for the Van der Heijden model. Takes two parameters h and
    k_c, with the latter being model-specific
    Returns the tuple (fn,ft) = (normal force, tangential force).
    """
    x = X[0]
    y = X[2]
    r_H = (x*x + y*y)**0.5
    delta_H = (r_H-h)*(tanh(1000*(r_H-h))+1)/2
    return (k_c*delta_H,0)
