import numpy as np

def buoy_dgl(dgl,eos,alpha = None,sbeta = None, g = None, lat = None,
        Tref = None, Sref = None):
    """Calculate the buoyancy field for a structured array dgl"""
    import mit
    b_exists = mit.check_field(dgl,'b')
    if not b_exists:
        dgl = mit.add_field(dgl,'b')
    S_exists = mit.check_field(dgl,'s')
    T_exists = mit.check_field(dgl,'t')
    if eos == 'lin':
        alpha, sbeta, g, Tref, Sref = lin_eos_params(alpha,sbeta,g,lat,Tref,Sref)
        for rw in np.arange(0,len(dgl[0]['y'])):
            if T_exists and not S_exists:
                """Thermal contribution only"""
                dgl[0]['b'][rw,:,:,:] = g*alpha*(dgl[0]['t'][rw,:,:,:] - Tref)
            if S_exists and not T_exists:
                """Haline contribution only"""
                dgl[0]['b'][rw,:,:,:] = -g*sbeta*(dgl[0]['s'][rw,:,:,:] - Sref)
            if T_exists and S_exists:
                dgl[0]['b'][rw,:,:,:] = g*alpha*(dgl[0]['t'][rw,:,:,:] - Tref)
                -g*sbeta*(dgl[0]['s'][rw,:,:,:] - Sref)

    elif eos == 'jmd':
        g = 9.81
        import jmd95
        import seawater
        import mit
        if not T_exists:
            T = np.zeros(np.shape(dgl[0]['s']))
        if not S_exists:
            S = np.zeros(np.shape(dgl[0]['t']))
        rhoConst = 1000
        ylen, xlen, zlen, tlen = mit.dgl_dims(dgl)
        #Pressure needed in decibars
        press = seawater.pres(-dgl[0]['z'],lat = 45)
        buoy_fac = g/rhoConst
        #Loop this routine along t and y to avoid running out of memory
        for ti in np.arange(0,tlen):
            for rw in np.arange(0,ylen):
                dgl[0]['b'][rw,:,:,ti] =-buoy_fac*(jmd95.densjmd95(
                np.squeeze(dgl[0]['s'][rw,:,:,ti]),
                np.squeeze(dgl[0]['t'][rw,:,:,ti]),press) - 1030)
    return dgl




def lin_eos_params(alpha = None,sbeta = None, g = None, lat = None,
        Tref = None, Sref = None):

        if g is None:
            g = 9.81
        if alpha is None:
            alpha = 2e-4
        if sbeta is None:
            sbeta = 7.4e-4
        if Tref is None:
            Tref = 0
        if Sref is None:
            Sref = 0
        return alpha, sbeta, g, Tref, Sref


def buoy_calc(Z, EOS, T = None, S = None, alpha = None,
              sbeta = None, g = None, lat = None,
              Tref = None, Sref = None):
    """Calculate the buoyancy given the temperature T, salinity S,depth,
    specified equation of state, thermal expansion coefficient (if necessary)
    haline coefficient (if necessary) and latitude (if necessary)"""
    if g is None:
        g = 9.81
    if EOS == 'lin':
        """Set default parameters if not otherwise specified"""
        if g is None:
            g = 9.81
        if alpha is None:
            alpha = 2e-4
        if sbeta is None:
            sbeta = 7.4e-4
        if Tref is None:
            Tref = 0
        if Sref is None:
            Sref = 0
        #Calculate the buoyancy contributions where relevant
        if S is None and T is not None:
            """Thermal contribution only"""
            buoy = g*alpha*(T - Tref)
        if T is None and S is not None:
            """Haline contribution only"""
            buoy = -g*sbeta*(S - Sref)
        if T is not None and S is not None:
            buoy = g*alpha*(T - Tref) -g*sbeta*(S - Sref)

    if EOS == 'jmd':
        import jmd95
        import seawater
        if T is None:
            T = np.zeros(np.shape(S))
        if S is None:
            S = np.zeros(np.shape(T))
        rhoConst = 1000
        #Pressure needed in decibars
        press = seawater.pres(-Z,45)
        #Loop this routine along t and y to avoid running out of memory
        sp = np.shape(S)
        buoy = np.empty(np.shape(S))
        for ti in np.arange(0,sp[-1]):
            for rw in np.arange(0,sp[0]):
                buoy_fac = g/rhoConst
                b =-buoy_fac*(jmd95.densjmd95(np.squeeze(S[rw,:,:,ti]),
                np.squeeze(T[rw,:,:,ti]),press) - 1030)
                buoy[rw,:,:,ti] = b
    return buoy
