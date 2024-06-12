#The Simulation of Optical Properties of Multilayer Structures by the Admittance Method
# definition of the function which calculates the coefficient R and T .
def coefficient_I(struct, wavelength, incidence, polarization):
    #n,d,lam,theta0):
    """
    This function computes the reflection and transmission coefficients
    of the structure using the fast impedance formalism.

    Args:
        struct (Structure): belongs to the Structure class
        wavelength (float): wavelength of the incidence light (in nm)
        incidence (float): incidence angle in radians
        polarization (float): 0 for TE, 1 (or anything) for TM

    returns:
        r (complex): reflection coefficient, phase origin at first interface
        t (complex): transmission coefficient
        R (float): Reflectance (energy reflection)
        T (float): Transmittance (energie transmission)


    R and T are the energy coefficients (real quantities)

    .. warning: The transmission coefficients have a meaning only if the lower medium
    is lossless, or they have no true meaning.
    """

    # In order to get a phase that corresponds to the expected reflected coefficient,
    # we make the height of the upper (lossless) medium vanish. It changes only the
    # phase of the reflection coefficient.

    # The medium may be dispersive. The permittivity and permability of each
    # layer has to be computed each time.
    if (struct.unit != "nm"):
        wavelength = conv_to_nm(wavelength, struct.unit)
    Epsilon, Mu = struct.polarizability(wavelength)
    thickness = np.array(copy.deepcopy(struct.thickness))
    # In order to ensure that the phase reference is at the beginning
    # of the first layer.
    thickness[0] = 0# wavelength/100
    Type = struct.layer_type
    # The boundary conditions will change when the polarization changes.
    # Wavevector in vacuum.
    k0 = 2 * np.pi / wavelength
    # Number of layers
    g = len(struct.layer_type)
    
    # Wavevector k_x, horizontal
    alpha = np.sqrt(Epsilon[Type[0]] * Mu[Type[0]]) * k0 * np.sin(incidence)
    # Computation of the vertical wavevectors k_z
    gamma = np.array(np.sqrt(
        Epsilon[Type] * Mu[Type] * k0 ** 2 - np.ones(g) * alpha ** 2))
    
    # Be cautious if the upper medium is a negative index one.
    if np.real(Epsilon[Type[0]]) < 0 and np.real(Mu[Type[0]]) < 0:
        gamma[0] = -gamma[0]

    # Changing the determination of the square root to achieve perfect stability
    if g > 2:
        gamma[1:g - 2] = gamma[1:g - 2] * (
                    1 - 2 * (np.imag(gamma[1:g - 2]) < 0))
    # Outgoing wave condition for the last medium
    if np.real(Epsilon[Type[g - 1]]) < 0 and np.real(
        Mu[Type[g - 1]]) < 0 and np.real(np.sqrt(Epsilon[
             Type[g - 1]] * Mu[Type[g - 1]] * k0 ** 2 - alpha ** 2)) != 0:
        gamma[g - 1] = -np.sqrt(Epsilon[Type[g - 1]] * Mu[
             Type[g - 1]] * k0 ** 2 - alpha ** 2)
    else:
        gamma[g - 1] = np.sqrt(Epsilon[Type[g - 1]] * Mu[Type[g - 1]] * k0 ** 2 - alpha ** 2)
        
    #to avoid divergence for large angles 
    if incidence>np.pi / 2 :
       incidence= -incidence+np.pi
    else:
       incidence=incidence
        
    # Effective reflective index TE

    n_s = np.zeros(g,dtype=complex)
    n_p = np.zeros(g,dtype=complex)
    n = np.array(np.sqrt(Epsilon*Mu))
    
    # Effective reflective index TE
    
    n_s[0] = n[Type[0]] * np.cos(incidence)
    n_s[1:] = np.sqrt(n[Type[1:]]**2 - n[Type[0]]**2 * np.sin(incidence)**2)
    opp = np.imag(n_s) > 0
    n_s = n_s - 2* n_s * (opp)
    
     # Effective reflective index TM
    
    n_p = n[Type]**2 / n_s
    
    # calculation of delta 
    
    delta = np.array(2*np.pi*thickness*n_s/wavelength)
    temp = -1.j*np.tan(delta)
    
    # choice of polarization for Initienization of the admittance  
    
    if polarization == 0:
        admittance = n_s
        Y = admittance[-1]
    else:
        admittance = n_p
        Y = admittance[-1]
    
    # Calculation of the complex admittence 
    
    PR =1
    for m in np.arange(g-2,-1,-1):    
     Y = (Y + admittance[m+1]*temp[m+1])/(1 + Y*temp[m+1]/admittance[m+1])
     if (m != 0):
    # calculation of product    
           PR *= (np.cos(delta[m]) - 1.j * Y * np.sin(delta[m])/admittance[m])

     # Calculation of r and t coefficient
    r = (admittance[0]-Y) / (admittance[0]+Y)
    t = (r+1)/ (PR)
    if (polarization == 1):
      r=-r
    # Calculation of R and T coefficient 
    R = abs(r)**2
    if polarization==0:
     admittance = n_s
     T = (admittance[-1].real/ admittance[0].real) * abs(t)**2  
     
    else:
      admittance = n_p
      T = (admittance[-1].real/ admittance[0].real) * abs(t)**2
    return(r,t,R ,T)
