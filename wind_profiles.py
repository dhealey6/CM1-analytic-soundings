import numpy as np
import math



def str_wind_profile(z,u_sfc,v_sfc,one_km_shear,three_km_shear,six_km_shear,dz,one_km_only,one_and_three_km,three_km_only,
                    one_three_six_km,six_km_only,three_and_six_km):
    #Inputs:
    #z = height array
    #u_sfc = zonal wind in m/s
    #v_sfc = meridional wind in m/s
    #one_km_shear = duh
    #three_km_shear = duh
    #six_km_shear = duh
    #dz = height spacing
    #one_km_only = user specified variable calling for only one shear layer in the lowest 1km.
    #one_and_three_km = user specified variable calling for two shear layers up to 1 and 3km AGL. 
    #three_km_only = as above, but calling for only one shear layer in the lowest 3km.
    #one_three_six_km = as above, but calling for three shear layers up to 1, 3, and 6 km AGL. 
    #six_km_only = as above, but calling for only one shear layer in the lowest 6km. 
    #three_and_six_km = as above, but calling for two shear layers up to 3 and 6km AGL. 
    
    #Create empty u and v wind arrays to be filled. These will be returned from the code.
    uwind = []    
    vwind = []
    
    #Find the indices for the specified heights here: 
    half_km_height = (np.where(z <= 500))[-1][-1]
    one_km_height = (np.where(z <= 1000))[-1][-1]
    two_km_height = (np.where(z <= 2000))[-1][-1]
    three_km_height = (np.where(z <= 3000))[-1][-1]
    five_km_height = (np.where(z <= 5000))[-1][-1]
    six_km_height = (np.where(z <= 6000))[-1][-1]
    eight_km_height = (np.where(z <= 8000))[-1][-1]
    
    #Add surface winds to the wind arrays first. 
    uwind.append(u_sfc)
    vwind.append(v_sfc)
    
    #Define actual shear value in s^-1. 
    shear1km = one_km_shear/1000
    shear3km = three_km_shear/3000
    shear6km = six_km_shear/6000
    
    #For only one shear layer in the lowest km. 
    if one_km_only:
        #looping through the lowest 1km:
        for h in z[:one_km_height]:
            uwind.append(uwind[-1] + ((shear1km)*(dz)))  #add u component based on shear specified. 
            vwind.append(vwind[-1])   #straight hodograph so v-component kept constant. 
        for h in z[one_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])            
    
    #For two shear layers, from 0-1 km AGL, and from 0-3 km AGL, keeping the 0-1 km layer as specified. 
    if one_and_three_km:
        #looping through lowest 1km:
        for h in z[:one_km_height]:
            uwind.append(uwind[-1] + ((shear1km)*(dz))) 
            vwind.append(vwind[-1])
        #looping through heights above 1km and below 3km:
        for h in z[one_km_height:three_km_height]:
            uwind.append(uwind[-1] + (((three_km_shear - one_km_shear)/(3000 - 1000))*(dz)))
            vwind.append(vwind[-1])
        #above three km: keep constant. 
        for h in z[three_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])
            
    #For only one shear layer in the lowest 3km.         
    if three_km_only:
        for h in z[:three_km_height]:
            uwind.append(uwind[-1] + ((shear3km)*(dz)))
            vwind.append(vwind[-1])
        for h in z[three_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])
            
    #For three shear layers, 0-1 km, 0-3 km, 0-6 km, keeping all values as specified. 
    if one_three_six_km:
        for h in z[:one_km_height]:
            uwind.append(uwind[-1] + ((shear1km)*(dz))) 
            vwind.append(vwind[-1])
        for h in z[one_km_height:three_km_height]:
            uwind.append(uwind[-1] + (((three_km_shear - one_km_shear)/(3000 - 1000))*(dz)))
            vwind.append(vwind[-1]) 
        for h in z[three_km_height:six_km_height]:
            uwind.append(uwind[-1] + (((six_km_shear - three_km_shear)/(6000 - 3000))*(dz)))
            vwind.append(vwind[-1])             
        for h in z[six_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])

    #For only one shear layer in the lowest 6km.         
    if six_km_only:
        for h in z[:six_km_height]:
            uwind.append(uwind[-1] + ((shear6km)*(dz)))
            vwind.append(vwind[-1])
        for h in z[six_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])        
            
    #For two shear layers, 0-3 km, 0-6 km, keeping all values as specified.         
    if three_and_six_km:
        for h in z[:three_km_height]:
            uwind.append(uwind[-1] + ((shear3km)*(dz))) 
            vwind.append(vwind[-1])
        for h in z[three_km_height:six_km_height]:
            uwind.append(uwind[-1] + (((six_km_shear - three_km_shear)/(6000 - 3000))*(dz)))
            vwind.append(vwind[-1])    
        for h in z[six_km_height+1:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])            
                    
    return uwind,vwind


def qc_wind_profile(qc_one_km_shear,qc_two_km_shear,qc_three_km_shear,height_curv,z,u_sfc,v_sfc,angle,dz):
    
    #Inputs:
    #qc_one_km_shear = one_km_shear specified. 
    #qc_two_km_shear = two_km_shear specified. 
    #qc_three_km_shear = three_km_shear specified. 
    #height_curv = height to which the quarter-circle will curve to. 
    #z = height array. 
    #u_sfc = zonal surface wind (in m/s)
    #v_sfc = meridional surface wind (in m/s)
    #angle = parameter for setting the size/radius of the quarter-circle portion of the hodograph. 
    #dz = vertical spacing. 
    
    #Create empty u and v wind arrays to be filled. These will be returned from the code. 
    uwind = []
    vwind = []
        
    #Find the indices for the specified heights here:     
    half_km_height = (np.where(z <= 500))[-1][-1]
    one_km_height = (np.where(z <= 1000))[-1][-1]
    two_km_height = (np.where(z <= 2000))[-1][-1]
    three_km_height = (np.where(z <= 3000))[-1][-1]
    five_km_height = (np.where(z <= 5000))[-1][-1]
    six_km_height = (np.where(z <= 6000))[-1][-1]
    eight_km_height = (np.where(z <= 8000))[-1][-1] 

    curvature_height = (np.where(z <= height_curv))[-1][-1]    

    #This is getting the indices in the height array to perform the curvature over.
    height_arr_indices = np.where(z <= height_curv)
    my_radians = np.linspace(0,(math.pi)/2,height_arr_indices[-1][-1]+1)
    x = 0
    
    #Getting the angle into the right units for the numpy math.
    new_angle = angle*((math.pi)/180)
    
    #Finding the u and v components of the one and two km shear vectors. These will be used when computing u and v components of the 
    #winds above the curvature height. 
    one_km_u = (qc_one_km_shear)*(np.cos(new_angle))
    one_km_v = (qc_one_km_shear)*(np.sin(new_angle))
    two_km_u = (qc_two_km_shear)*(np.cos(new_angle))
    two_km_v = (qc_two_km_shear)*(np.sin(new_angle))    
    
    #Setting the curvature shear vector components equal to the v component of whichever height the curvature will remain beneath. 
    #This keeps the curvature in a nice quarter-circle. 
    curvature_u_one = one_km_v
    curvature_v_one = one_km_v
    curvature_u_two = two_km_v
    curvature_v_two = two_km_v
    
    #Finding the u component of the three km shear vector for calculating the wind above the curvature height up to 3km. 
    three_km_shear_u = ((qc_three_km_shear**2) - (one_km_v**2))**(1/2)
     
    #For curvatures below 1km, use this loop:     
    if height_curv < 1000:
        for h in z[:curvature_height]:
            uwind.append(curvature_u_one*(1-np.cos(my_radians[x]))+u_sfc)
            vwind.append(curvature_u_one*(np.sin(my_radians[x]))+v_sfc)            
            x = x+1            
        for h in z[curvature_height:one_km_height]:
            uwind.append(uwind[-1] + ((one_km_u - curvature_u_one)/(1000 - height_curv))*(dz))
            vwind.append(vwind[-1] + ((one_km_v - curvature_v_one)/(1000 - height_curv))*(dz))
        for h in z[one_km_height:three_km_height]:
            uwind.append(uwind[-1] + ((three_km_shear_u - one_km_u)/(3000 - 1000))*(dz))
            vwind.append(vwind[-1])
        for h in z[three_km_height:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])
    
    #For curvatures above or equal to 1km, and less than 2km, use this loop: 
    if height_curv >= 1000 and height_curv < 2000:
        for h in z[:curvature_height]:
            uwind.append(curvature_u_two*(1-np.cos(my_radians[x]))+u_sfc)
            vwind.append(curvature_u_two*(np.sin(my_radians[x]))+v_sfc)            
            x = x+1            
        for h in z[curvature_height:two_km_height]:
            uwind.append(uwind[-1] + ((two_km_u - curvature_u_two)/(2000 - height_curv))*(dz))
            vwind.append(vwind[-1] + ((two_km_v - curvature_v_two)/(2000 - height_curv))*(dz))
        for h in z[two_km_height:three_km_height]:
            uwind.append(uwind[-1] + ((three_km_shear_u - two_km_u)/(3000 - 2000))*(dz))
            vwind.append(vwind[-1])
        for h in z[three_km_height:]:
            uwind.append(uwind[-1])
            vwind.append(vwind[-1])            

    
    return uwind, vwind