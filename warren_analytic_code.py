import numpy as np
from metpy.calc import (saturation_vapor_pressure,exner_function,vapor_pressure,dewpoint,
                       relative_humidity_from_dewpoint,saturation_mixing_ratio,potential_temperature,
                       wind_direction,wind_speed,bunkers_storm_motion)
from metpy.units import units
import math

import sharppy
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params


def warren_analytic_sounding(z_mlt,p_mlt,t_mlt,rh_mlt,lrth_ml,cape0,z_bmax,z_tp,rh_tp,rhfact,g,p00,cpd,cpv,cpl,rd,rv,
                             eps,t0,lv0,lv1,lv2,capetol,ztop,dz,nz,z,wind_spd,wind_dir):

    #Create arrays for environmental profiles:
    t = np.zeros(nz)     #temperature
    td = np.zeros(nz)    #dewpoint
    tv = np.zeros(nz)    #virtual temperature
    th = np.zeros(nz)    #potential temperature
    thv = np.zeros(nz)   #virtual potential temperature
    qv = np.zeros(nz)    #mixing ratio
    pi = np.zeros(nz)    #Exner function
    p = np.zeros(nz)     #pressure
    rh = np.zeros(nz)    #relative humidity

    #Create arrays for parcel profiles:
    t_p = np.zeros(nz)   #temperature of parcel
    qv_p = np.zeros(nz)  #mixing ratio of parcel
    tv_p = np.zeros(nz)  #virtual temperature of parcel

    #Starting to create the profile: McCaul and Weisman (2001) approach with the code from Rob Warren
    #Warren et al. (2017) 

    #MIXED-LAYER:

    #Find the top index of the mixed layer: 
    ind1 = np.where(z <= z_mlt)
    k1 = ind1[-1][-1] 

    #store variables at mixed layer top:
    t[k1] = t_mlt
    p[k1] = p_mlt

    #determine mixing ratio at mixed layer top:
    e = rh_mlt*(((saturation_vapor_pressure(t_mlt*units.kelvin)).magnitude)*100)
    qv[k1] = (eps*e)/(p[k1]-e)

    #determine nondimensional pressure (Exner function) at mixed layer top:
    pi[k1] = exner_function(p[k1]*units.Pa)

    #determine potential temperature at mixed layer top: 
    th[k1] = t[k1]/pi[k1]

    #determine the virtual potential temperature at mixed layer top: 
    thv[k1] = th[k1]*(1+qv[k1]/eps)/(1+qv[k1])

    #loop downward from mixed layer top:
    for k in range(k1-1,-1,-1):
        #update potential temperature:
        th[k] = th[k+1]-(lrth_ml/1000)*dz

        #set initial mixing ratio to that at the previous level:
        qv[k] = qv[k+1]

        #iterate to get the pressure, mixing ratio, and virtual temperature:
        for i in range(0,10):
            #qv0 = qv(k)

            #compute the virtual potential temperature:
            thv[k] = th[k]*(1+qv[k]/eps)/(1+qv[k])

            #compute the pressure:
            pi[k] = pi[k+1]+g*dz/(cpd*0.5*(thv[k]+thv[k+1]))
            p[k] = p00*(pi[k]**(cpd/rd))

            #compute the temperature:
            t[k] = th[k]*pi[k]

            #note the layer-mean temperature and mixing ratio:
            tbar = 0.5*(t[k]+t[k+1])
            qvbar = 0.5*(qv[k]+qv[k+1])

            #compute the layer-mean latent heat:
            lv = lv1 - lv2*tbar

            #compute the gas constant and specific heat for the layer:
            rm = rd+rv*qvbar
            cpm = cpd+cpv*qvbar

            #update the mixing ratio:
            qv[k] = qv[k+1]+cpm*tbar/lv*(rm/cpm*math.log(p[k]/p[k+1])-math.log(t[k]/t[k+1]))

    #compute the nondimensional pressure (exner function) profile:
    pi[ind1] = exner_function(p[ind1]*units.Pa)

    #compute the potential temperature profile:
    th[ind1] = t[ind1]/pi[ind1]

    #compute the virtual temperature profile:
    tv[ind1] = t[ind1]*(1+qv[ind1]/eps)/(1+qv[ind1])

    #compute the dewpoint temperature profile:
    td[ind1] = (dewpoint(vapor_pressure(p[ind1]*units.Pa,qv[ind1])).magnitude)+273.15

    #compute the relative humidity profile:
    rh[ind1] = relative_humidity_from_dewpoint(t[ind1]*units.kelvin,td[ind1]*units.kelvin).magnitude

    #set parcel properties below the mixed layer top: 
    t_p[ind1] = th[k1]*pi[ind1]
    qv_p[ind1] = saturation_mixing_ratio(p_mlt*units.Pa,t_mlt*units.kelvin)
    tv_p[ind1] = t_p[ind1]*(1+qv_p[ind1]/eps)/(1+qv_p[ind1])
    
    #MIXED LAYER TOP TO TROPOPAUSE

    #Find the levels between mixed layer top and tropopause:
    ind2 = np.where((z >= z_mlt) & (z < z_tp))
    k2 = ind2[-1][-1]

    #compute relative humidity profile:
    rh[ind2] = rh_mlt-(rh_mlt-rh_tp)*((z[ind2]-z_mlt)/(z_tp-z_mlt))**rhfact

    #if I wanted to add a dry layer, here is where I would do it. But for now, will not include.

    #compute buoyancy profile:
    zprime = z[ind2]-z_mlt
    zprime_bmax = z_bmax-z_mlt
    b = cape0*zprime/(zprime_bmax**2)*np.exp(-1*zprime/zprime_bmax)

    kk3 = len(ind2[-1])
    kk4 = kk3-2
    kk5 = kk3-1

    #compute the mixed-layer top to tropopause CAPE:
    cape = np.sum(0.5*(b[0:kk4]+b[1:kk5])*dz)

    #add in a restart?

    #scale buoyancy to account for the difference:
    b = b*cape0/cape

    #loop over levels:
    for k in range(k1+1,k2+1):
        #set initial virtual temperature to that at previous level:
        tv[k] = tv[k-1]

        #set the initial parcel values for the iteration:
        th1 = potential_temperature(p[k-1]*units.Pa,t_p[k-1]*units.kelvin).magnitude
        t1 = t_p[k-1]
        qv1 = qv_p[k-1]
        thlast=th1

        #iterate to get the pressure, virtual temperature, and parcel properties:
        converged=0
        count=0
        while converged == 0:
            #compute the pressure:
            p[k] = p[k-1]*math.exp(-1*g*dz/(rd*0.5*(tv[k-1]+tv[k])))

            #compute the parcel temperature:
            t2 = thlast*exner_function(p[k]*units.Pa)

            #compute the parcel mixing ratio:
            qv2 = saturation_mixing_ratio(p[k]*units.Pa,t2*units.kelvin)

            #compute the mean parcel temperature and mixing ratio for the layer:
            tbar = 0.5*(t1+t2)
            qvbar = 0.5*(qv1+qv2)

            #compute the latent heat of vaporization for the layer:
            lv = lv1-lv2*tbar

            #compute the gas constant and total specific heat for the layer:
            rm = rd+rv*qvbar
            cpm = cpd+cpv*qvbar

            #update the parcel potential temperature:
            th2 = th1*math.exp((rm/cpm-rd/cpd)*math.log(p[k]/p[k-1])-lv/(cpm*tbar)*(qv2-qv1))

            #check for convergence:
            if abs(th2-thlast) >= 0.0002:
                thlast = thlast+0.3*(th2-thlast)
            else:
                converged = 1

            #compute the parcel virtual temperature:
            tv2 = t2*(1+qv2/eps)/(1+qv2)

            #compute the environment virtual temperature:
            tv[k] = g*tv2/(b[k-k1]+g)

            count = count+1

        #store parcel properties:
        t_p[k] = t2
        qv_p[k] = qv2
        tv_p[k] = tv2

        #set the initial temperature to the virtual temperature:
        t[k] = tv[k]

        #iterate to get the temperature and mixing ratio:
        for i in range(1,20):
            #compute the mixing ratio:
            e = rh[k]*(((saturation_vapor_pressure(t[k]*units.kelvin)).magnitude)*100)
            qv[k] = eps*e/(p[k]-e)

            #compute the temperature:
            t[k] = tv[k]*(1+qv[k])/(1+qv[k]/eps)        

    qv[k2] = qv[k2-1]        

    #compute the exner pressure:
    pi[ind2] = exner_function(p[ind2]*units.Pa)

    #compute theta:
    th[ind2] = potential_temperature(p[ind2]*units.Pa,t[ind2]*units.kelvin).magnitude

    #compute theta-v:
    thv[ind2] = th[ind2]*(1+qv[ind2]/eps)/(1+qv[ind2])

    #compute dewpoint profile:
    td[ind2] = (dewpoint(vapor_pressure(p[ind2]*units.Pa,qv[ind2])).magnitude)+273.15

    #TROPOPAUSE TO MODEL TOP

    #Find levels above the tropopause and below model top:
    ind3 = np.where(z >= z_tp)
    k3 = ind3[-1][-1]

    #set temperature constant:
    t[ind3] = t[k2]

    #set mixing ratio constant:
    qv[ind3] = qv[k2]

    #compute virtual temperature:
    tv[ind3] = t[ind3]*(1+qv[ind3]/eps)/(1+qv[ind3])

    #compute the pressure profile:
    for k in range(k2+1,k3+1):
        p[k] = p[k-1]*np.exp(-1*g*dz/(rd*0.5*(tv[k-1]+tv[k])))

    #get dewpoint profile:
    #td[ind3] = (dewpoint(vapor_pressure(p[ind3]*units.Pa,qv[ind3])).magnitude)+273.15
    td[ind3] = td[k2-1]

    #add the iteration for cape here:


    #compute the nondimensional pressure (exner function):
    pi[ind3] = exner_function(p[ind3]*units.Pa)

    #compute potential temperature:
    th[ind3] = potential_temperature(p[ind3]*units.Pa,t[ind3]*units.kelvin).magnitude

    #compute theta-v:
    thv[ind3] = th[ind3]*(1+qv[ind3]/eps)/(1+qv[ind3])

    #get relative humidity profile:
    rh[ind3] = relative_humidity_from_dewpoint(t[ind3]*units.kelvin,td[ind3]*units.kelvin).magnitude

    #set parcel properties above tropopause: 
    t_p[ind3] = (potential_temperature(p[k2]*units.Pa,t_p[k2]*units.kelvin).magnitude)*pi[ind3]
    qv_p[ind3] = saturation_mixing_ratio(p[k2]*units.Pa,t[k2]*units.kelvin)
    tv_p[ind3] = t_p[ind3]*(1+qv_p[ind3]/eps)/(1+qv_p[ind3])

    #compute maximum lapse rate:
    dthdz = (th[1:nz-1]-th[0:nz-2])/dz

    #check for superadiabatic lapse rates:
    if min(dthdz) <= 0:
        print("PROFILE CONTAINS SUPERADIABATIC LAPSE RATE!")
        
    return t, td, tv, th, thv, qv, pi, p, rh, t_p, qv_p, tv_p
    