# -*- coding: utf-8 -*-
"""

Author: Zeke Noller

READ ME: You can edit the initial parameters between the asterisks. You can write to\
the filepath by changing the variable 'write_to_excel' to True


"""

import numpy as np
import math 
import scipy
import scipy.sparse
import matplotlib.pyplot as plt 

# **************************************************

## Initial structural parameters
t_skin = 5e-3 #m, Skin thickness
t_web = 3e-3 #m, Web thickness
t_sparcap = 3e-3 #m, Sparcap material thickness
t_str = 3e-3 #m, Stringer material thickness
l_sparcap = 25e-3 #m, Sparcap leg length 
l_str = 25e-3 #m, Stringer leg length 

lam = 10

write_to_excel = False

# **************************************************

try:
    import IPython
    shell = IPython.get_ipython()
    shell.enable_matplotlib(gui='inline') #qt5, inline
except:
    print('Unable to open plotting window')

## Airfoil shape
import geometry 

### Description of system architecture and algorithm ###
## Inputs: 
# x = vector of decision variables
# f = vector of function variables 
# g = vector of constraint variables 

# i = counting variable for x variables
# j = counting variable for constraint functions (g)

# m = number of x variables 
# n = number of constraint functions (g)

## Conditions: 
# len(x) > len(g) i.e. there are more variables than constraint functions

## Assumptions (for aerofoil problem)
# f is function to be minimised: cross-sectional area 
# len(f) = 1
# f has minimum and is relatively convex 
# g contains stresses - bending, shear at each critical element 
# g goes to infinity at x_vec = 0 vector

## Material properties
## Al-7075-T73 from MIL handbook (2003) 3.7.6.0(b3) pp.3-373
## Sheet S-basis ##
Ftu_sh = 67*6.89e6 #Pa, ultimate tensile
Fty_sh = 56*6.89e6 #Pa, yield tensile
Fcy_sh = 55*6.89e6 #Pa, yield compressive
Fsu_sh = 38*6.89e6 #Pa, shear ultimate
E_sh = 10.3*6.89e6 #Pa, tensile modulus
Ec_sh = 10.5*6.89e9 #Pa, compressive modulus
Gsh = 3.9*6.89e9 #Pa, shear modulus
miu_sh = 0.33 # Poisson's ratio (dimensionless)

## Al-7075-T73 from MIL handbook 3.7.6.0(g2) pp.3-3386
## Extrusion A-basis ##
Ftu_ex = 66*6.89e6 #Pa, ultimate tensile
Fty_ex = 56*6.89e6 #Pa, yield tensile
Fcy_ex = 58*6.89e6 #Pa, yield compressive
Fsu_ex = 37*6.89e6 #Pa, shear ultimate
E_ex = 10.4*6.89e6 #Pa, tensile modulus
Ec_ex = 10.7*6.89e9 #Pa, compressive modulus
Gex = 4.0*6.89e9 #Pa, shear modulus
miu_ex = 0.33 # Poisson's ratio (dimensionless)

## Wing geometry parameters (constant)
span = geometry.span
semispan = geometry.semispan
planform_area = geometry.planform_area
chord = geometry.chord

## Absolute spar coordinates. Arbitrarily chosen
left_spar_x = 0.20889877 #m
upper_left_spar_z = 0.108301478 #m
lower_left_spar_z = -0.065517339 #m

right_spar_x = 1.211612868 #m
upper_right_spar_z = 0.089759624 #m
lower_right_spar_z = -0.038556028 #m

wingbox_length = right_spar_x - left_spar_x

## Heights of spars
right_dH = (upper_right_spar_z - lower_right_spar_z) # m
left_dH = (upper_left_spar_z - lower_left_spar_z) # m
sigma_dH = right_dH + left_dH

## For centroid calculations
Cx = (right_spar_x + left_spar_x)/2
Cz = (upper_left_spar_z + lower_left_spar_z + \
      upper_right_spar_z + lower_right_spar_z)/4

## For interpolation 
upper_mid_x = 0.626696311 #m
upper_mid_z = 0.134236679 #m

lower_mid_x = upper_mid_x
lower_mid_z = -0.085329299 #m

# Note these are clockwise from the origin, including the origin again.
OML_x = geometry.x_outer
OML_z = geometry.z_outer
OML_len = len(OML_x)
half_OML_len = OML_len/2

# ! What about the concavity of the aerofoil - the circumference will be 
# underestimated
S = 0
for count in range(OML_len):
    if count < OML_len - 1:
        S += np.sqrt((OML_x[count+1] - OML_x[count])**2 + (OML_z[count+1] - OML_z[count])**2)

# Split coordinates. Clockwise from top left
OML_x_upper = geometry.x_outer_upper
OML_x_lower = geometry.x_outer_lower
OML_z_upper = geometry.z_outer_upper
OML_z_lower = geometry.z_outer_lower

left_right_OML_z = geometry.left_right_OML_z
left_right_OML_z = geometry.left_right_OML_z


str_mult = 1.5
n_str = int(wingbox_length / (str_mult*l_str))

h = np.array([3e-5, # t_skin
              3e-5, # t_web
              5e-5, # t_sparcap
              5e-5, # t_str
              5e-5, # l_sparcap
              5e-5, # l_str
              2 # n_str - must be whole number to make sense. Equal on upper and lower surfaces 
              ])

h_mat = np.diag(h)

# init_vec = np.array((t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str,a_rib))
init_vec = np.array((t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str))
with np.printoptions(precision=3, suppress=True):
    print(f"The initial parameters are: {init_vec}")

m = len(init_vec) # Number of decision variables. Counting variable 'i'

n = 3 # Number of constraint functions. Counting variable 'j'

## Initialising linear algebra 
# In form A @ x_vec_long = b
A = np.zeros((m,m)) # Contains del2F and delGj
x_vec_long = np.zeros((m,1)) # Contains delta_x and delta_lambda i.e. solution step
b = np.zeros((m,1)) # Contains del1F and G vector

Z_stat = 2

wingbox_area = 0

## Initialising global parameters

Ce_1 = 0.317 # Crippling coefficient: aluminium, 1 free edge
Ce_2 = 0.295 # Crippling coefficient: aluminium, 2 free edge

## Forces, torques and moments 
# Pull-up manoeuvre
MTOW = 3620 #kg, maximum take-off weight
acceleration = 2.5 #g
gravity = 9.81 #m/s2
Fz = -MTOW*acceleration*gravity/2 # divide by two because two wings
Mxx = Fz*semispan # maximum bending moment 
T = -Fz * (Cx - 0.25*chord) # positive clockwise

## NOTE: signs of Fz and Mxx are by convention. In reality the loads are vertical

# ! Change to be boundary function i.e. f - lambda*log(G(x))
# ! Careful of too large a step to the left - can get into infeasible i.e. negative
# ! parameters really quickly 

## Updates grad vector of total area 
def update_del1F(x_vec):
    # ! Must be before delGj derivatives are calculated as they are for variables away from current x_vec
    # Hence unpacking variables from x_vec instead of taking global values for x_i

    t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str = x_vec
    
    global b # Tells Python to let update_del1F to write changes to b
    b[0,0] = S # dF/d skin thickness
    b[1,0] = sigma_dH # dF/d web thickness
    b[2,0] = 8*l_sparcap - 8*t_sparcap # dF/d sparcap thickness
    b[3,0] = 3*n_str*l_str - 4*n_str*t_str # dF/d stringer thickness
    b[4,0] = 8*t_sparcap # dF/d sparcap leg length 
    b[5,0] = 3*n_str*t_str # dF/d stringer leg length 
    b[6,0] = 3*l_str*t_str - 2*t_str**2 # dF/d number of stringers
    
## Updates Hessian submatrix for second partial derivatives of total area
def update_del2F(x_vec):
    # ! Must be before delGj derivatives are calculated as they are for variables away from current x_vec
    # Hence unpacking variables from x_vec instead of taking global values for x_i

    t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str = x_vec
    
    global A # Tells Python to let update_del2F to write changes to A
    ## Note A is symmetric because of commutativity of differentiation
    A[2,2] = -8
    A[3,3] = -4*n_str
    A[4,2] = A[2,4] = 8
    A[5,3] = A[3,5] = 3*n_str
    A[6,3] = A[3,6] = 3*l_str - 4*t_str
    A[5,6] = A[6,5] = 3*t_str

def wingbox_area(x_vec):
 
    t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str = x_vec
    n_str = int(wingbox_length / (str_mult*l_str))  

    wingbox_area = 0
    wingbox_area += S*t_skin
    wingbox_area += (right_dH + left_dH) * t_web
    wingbox_area += 4 * (2*l_sparcap*t_sparcap - t_sparcap**2)
    wingbox_area += 2*n_str * (3*l_str*t_str - 2*t_str**2)
    
    return wingbox_area

## Returns central, finite, numerical derivative del(G_j) for each of the n stress margins examined 
def G(x_vec):
    
    ## Unpacking vector into local variables 
    t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str = x_vec
    n_str = int(wingbox_length / (str_mult*l_str)) 

    
    # Stringer coordinates
    stringer_x_vec = np.linspace(left_spar_x,right_spar_x,num=int(n_str)+2)
    stringer_lower_vec = np.zeros_like(stringer_x_vec)
    stringer_upper_vec = np.zeros_like(stringer_x_vec)
    
    for index, x in enumerate(stringer_x_vec):
        stringer_upper_vec[index] = np.interp(x,OML_x_upper,OML_z_upper)
    for index, x in enumerate(stringer_x_vec[::-1]):       
        stringer_lower_vec[index] = np.interp(x,OML_x_lower[::-1],OML_z_lower[::-1])
    stringer_lower_vec = stringer_lower_vec
            
    ## Stringers - Z-section - two L-sections 
    b1_str = l_str - 0.5*t_str
    b2_str = 0.5*l_str - 0.5*t_str
    b_dash_str = (b1_str + b2_str)/(2*t_str)
    str_cripple_stress = Ce_1 * ((Fcy_ex * Ec_ex)**0.5)/(b_dash_str**0.75)
    w_eff = 1.7*t_skin * (Ec_sh / str_cripple_stress)**0.5
    
    A_str = 3*l_str*t_str - 2*t_str**2
    
    A_sk = t_skin*w_eff
    
    ## Sparcaps - L-section 
    b1_sp = l_sparcap - 0.5*t_sparcap
    b_dash_sp = (b1_sp)/(t_sparcap)
    sp_cripple_stress = Ce_2 * ((Fcy_ex * Ec_ex)**0.5)/(b_dash_sp**0.75)
    
    A_sparcap = 2*l_sparcap*t_sparcap - t_sparcap**2
    
    ## Constructing calculation table
    n_columns = 21 # just a guess
    n_elements = int(4+2*n_str)
    n_elements_half = int(0.5*n_elements)
    global table_calc
    table_calc = np.zeros((n_elements,n_columns))

    ## Note: clockwise, with element 0 being top left and web above it 
    
    # Type: 1 = stringer, 2 = spar car
     
    ## Column 0: Element number
    for row in range(n_elements): table_calc[row,0] = row
    
    ## Column 1: Element labels
    table_calc[:,1] = 1 # Stringers
    table_calc[0,1] = table_calc[n_str+1,1] = \
        table_calc[n_str+2,1] = table_calc[-1,1] = 2 # Spar caps
        
    ## Column 2: X
    table_calc[:n_str+2,2] = stringer_x_vec
    table_calc[n_str+2:,2] = stringer_x_vec[::-1]
    
    ## Column 3: Z
    table_calc[:n_elements_half,3] = stringer_upper_vec 
    table_calc[0,3] += - 0.5*l_sparcap 
    table_calc[n_str+1,3] += - 0.5*l_sparcap
    table_calc[1:n_str+1,3] += - 0.5*l_str
    
    table_calc[n_str+2:,3] = stringer_lower_vec
    table_calc[n_str+2,3] += 0.5*l_sparcap
    table_calc[n_elements-1,3] += 0.5*l_sparcap
    table_calc[n_str+2:n_elements-1,3] += 0.5*l_str
            
    h_sparcap_left = table_calc[0,3] - table_calc[-1,3]
    h_sparcap_right = table_calc[n_str+1,3] - table_calc[n_str+2,3]
  
    ## Column 4: x
    table_calc[:,4] = table_calc[:,2] - Cx

    ## Column 5: z
    table_calc[:,5] = table_calc[:,3] - Cz
    
    ## Column 6: A_effective
    for stringer in range(1,n_str+1):
        h_st = stringer_upper_vec[stringer] - stringer_lower_vec[stringer]
        h_sk = h_st + l_str
        A_str_eff = A_sk + A_str*(h_st/h_sk)**2
        table_calc[stringer,6]
        table_calc[stringer,6] = table_calc[stringer + n_str + 2,6] = A_str_eff
    
    table_calc[0,6] = table_calc[n_str+1,6] = \
        A_sparcap + (t_web * left_dH / 6) * (left_dH / h_sparcap_left)**2
        
    table_calc[n_str+2,6] = table_calc[n_elements-1,6] = \
        A_sparcap + (t_web * right_dH / 6) * (right_dH / h_sparcap_right)**2
    
    ## Column 7: Ae*x**2
    table_calc[:,7] = table_calc[:,6] * table_calc[:,4]**2
    
    ## Column 8: Ae*z**2
    table_calc[:,8] = table_calc[:,6] * table_calc[:,5]**2

    ## Column 9: Ae*x*z
    table_calc[:,9] = table_calc[:,6] * table_calc[:,4] * table_calc[:,5]
    
    global Ixx, Izz, Ixz
    Izz = np.sum(table_calc[:,7]) 
    Ixx = np.sum(table_calc[:,8])
    Ixz = np.sum(table_calc[:,9])
    
    
    # Note Sz is -ZZ N because it's pointing up 
    # Mxx is -XX N.m because I said so 
    # Torque is positive in the clockwise direction 
    
    ## Column 10: Bending stresses
    # (xi(IxxMz - IzyMx) - zi(IzzMx - IxzMz)) / (IxxIzz - Ixz^2)
    # NOTE: Mz = 0
    # => (-xiIxzMx - ziIzzMx) / (IxxIzz - Ixz^2)
    table_calc[:,10] = (table_calc[:,4]*-Ixz*Mxx + table_calc[:,5]*Izz*Mxx) / \
        (Ixx*Izz - Ixz**2)
        
    ## Column 11: dx 
    for row in range(n_elements-1):
        table_calc[row,11] = table_calc[row+1,4] - table_calc[row,4]
    table_calc[-1,11] = 0
    
    ## Column 12: dz
    for row in range(n_elements-1):
        table_calc[row,12] = table_calc[row+1,5] - table_calc[row,5]
    table_calc[-1,12] = left_dH    
    
    ## Column 13: Enclosed area using Green's theorem 
    # A = 0.5 * Sum(xi*zi+1 - xi+1*zi)
    for row in range(0,n_elements-1):
        table_calc[row,13] = 0.5*abs((table_calc[row,4] * table_calc[row+1,5] - \
            table_calc[row+1,4] * table_calc[row,5]))
    table_calc[-1,13] = 0.5*abs((table_calc[-1,4]*table_calc[0,5] - table_calc[0,4]*table_calc[-1,5]))
    
    global A_wingbox
    A_wingbox = np.sum(table_calc[:,13])
    
    ## Column 14: q'
    # q' = Ai(xi(IxxSz - IxzSz) + zi(IzzSz - IxzSx))/(IxxIzz - Ixz^2) + q'i-1
    # Sx = 0 
    # => q' = Ai(-xi(IxzSz) + zi(IyzSz))/(IxxIzz - Ixz^2) + q'i-1

    # table_calc[:,10] = (table_calc[:,4]*-Ixz*Mxx - table_calc[:,5]*Izz*Mxx) / \
        # (Ixx*Izz - Ixz**2)
        
    for row in range(1,n_elements): 
        table_calc[row,14] = table_calc[row,6]*(table_calc[row,4]*-Ixz*Fz + table_calc[row,5]*Izz*Fz) / \
            (Ixx*Izz - Ixz**2)
        table_calc[row,14] += table_calc[row-1,14]
    
    ## Column 15: q'A
    table_calc[:,15] = table_calc[:,14] * table_calc[:,13]
      
    global sum_qdash_Ae
    sum_qdash_Ae = np.sum(table_calc[:,15])
    
    q1 = (-sum_qdash_Ae / A_wingbox) + T/(2 * A_wingbox)
        
    ## Column 16: q
    table_calc[:,16] = table_calc[:,14] + q1
    
    ## Column 17: Shear stresses
    for row in range(-1,n_elements-1):
        table_calc[row,17] = table_calc[row,16]
        if table_calc[row+1,1] == 1:
            table_calc[row,17] /= t_str 
        elif table_calc[row+1,17] == 2:
            table_calc[row,1] /= t_sparcap

    ## Column 18: Tensile margins of safety 
    ## NOTE: Using ultimate load because tensile yield stress is more than 2/3 of tensile ultimate stress
    for row in range(n_elements):
        bending_stress = table_calc[row,10]
        if bending_stress > 0:
            table_calc[row,18] = Ftu_ex/(1.5 * bending_stress) - 1
    
    MS_tens = table_calc[:,18]
    MS_tens_nonzero = MS_tens[np.nonzero(MS_tens)]
    G1 = np.min(MS_tens_nonzero)
    G1_av = np.average(MS_tens_nonzero)
    G1_std = np.std(MS_tens_nonzero)
    # G1 = G1_av - Z_stat*G1_std

    ## Column 19: Compressive margins of safety 
    ## NOTE: Using ultimate load because failure by crippling is ultimate
    for row in range(n_elements):
        bending_stress = table_calc[row,10]
        if bending_stress < 0:
            bending_stress = abs(bending_stress)
            if table_calc[row,1] == 1:
                table_calc[row,19] = str_cripple_stress/(1.5 * bending_stress) - 1
            elif table_calc[row,1] == 2:
                table_calc[row,19] = sp_cripple_stress/(1.5 * bending_stress) - 1
                
    MS_comp = table_calc[:,19]
    MS_comp_nonzero = MS_comp[np.nonzero(MS_comp)]
    G2 = np.min(MS_comp_nonzero)
    G2_av = np.average(MS_comp_nonzero)
    G2_std = np.std(MS_comp_nonzero)
    # G2 = G2_av - Z_stat*G2_std
    
    ## Column 20: Shear margins of safety 
    ## NOTE: Using ultimate load because no shear yield given
    for row in range(n_elements):
        shear_stress = abs(table_calc[row,17])
        table_calc[row,20] = Fsu_sh/(1.5 * shear_stress) - 1

    MS_shear = table_calc[:,20]
    MS_shear_nonzero = MS_shear[np.nonzero(MS_shear)]
    G3 = np.min(MS_shear_nonzero)
    G3_av = np.average(MS_shear_nonzero)
    G3_std = np.std(MS_shear_nonzero)
    # G3 = G3_av - Z_stat*G3_std
    
    ## Additional tweaking functions e.g. keep variables positive, keep stringers from overlapping
    
    output_vec = np.zeros((3+m,1))
    output_vec[:3,0] = np.array([G1,G2,G3])
    output_vec
    output_vec[3:] = 1000*x_vec.reshape(m,1)
    output_vec[-1] /= 100000
    

    return output_vec
    
#  0      1        2       3        4     5      6    7  (for reference)
# t_skin,t_web,t_sparcap,t_str,l_sparcap,l_str,n_str,a_rib 

def update_matrices(x_vec):
    global A, b
    
    update_del2F(x_vec) # Updates upper left of A matrix to del2F
    update_del1F(x_vec) # Updates first m entries of b vector to del1F
    
    # G_mat = np.zeros((n,m))
    for i in range(m): # All elements in square A matrix
        G_right = G(x_vec + h_mat[:,i])
        G_left = G(x_vec - h_mat[:,i])
        G0 = G(x_vec)
        G1_leftright = (G_right - G_left ) / (2*h[i])
        b[i,0] += -(1/lam)*np.sum(G1_leftright / G0)
        for i2 in range(m):
            if i == i2:
                G2_same = (G_left - 2*G0 + G_right) / (h[i]**2)
                A[i,i2] += (1/lam)*np.sum((G2_same/G0) - ((G1_leftright)**2)/(G0**2))
            else:
                G_up = G(x_vec + h_mat[:,i2])
                G_down = G(x_vec + h_mat[:,i2])
                G1_updown = (G_up - G_down) / (2*h[i2])
                G_topright = G(x_vec + h_mat[:,i] + h_mat[:,i2])
                G_bottomright = G(x_vec + h_mat[:,i] - h_mat[:,i2])
                G_bottomleft = G(x_vec - h_mat[:,i] - h_mat[:,i2])
                G_topleft = G(x_vec - h_mat[:,i] + h_mat[:,i2])
                
                G2_mixed = (G_topright + G_bottomleft - G_bottomright - G_topleft)/ \
                    (4*h[i]*h[i2])
                A[i,i2] += (1/lam)*np.sum((G2_mixed/G0) - (G1_leftright*G1_updown)/(G0**2))

x = np.ones((m,1))

num_iter = 16
twenty_iter = np.zeros((3+m+1,num_iter))
twenty_axis = np.linspace(0,num_iter-1,num=num_iter)

def iterate():

    global x, step
    
    x[:,0] = init_vec
    x_new = np.copy(x)
    step = np.ones_like(x)
    tol = 1e-6

    it = 0 

    while np.min(np.abs(step)) > tol and it < num_iter:

        twenty_iter[:-1,it] = G(x)[:,0]
        twenty_iter[-1,it] = 1000*wingbox_area(x) # For graphing multiply by 100
        print(it)
        
        if np.min(G(x+(1/lam)*step)) < 0 or np.min(x + (1/lam)*step) < 0:
            print(f"*"*50)
            print(f"The final solution is {x.T}")
            print(f"The final minimum factors of safety are {G(x).T}")
            print(f"The final area is {wingbox_area(x)}m^2")
            break

        if it == 0:
            update_matrices(init_vec)

            step = np.linalg.solve(A,-b)

            x += (1/lam)*step
            
            x[6,0] = int(x[6,0])


        else:
            update_matrices(x[:,0])

            step = np.linalg.solve(A,-b)

            x += (1/lam)*step

            x[6,0] = int(wingbox_length / (str_mult*x[4,0]))
            
        with np.printoptions(precision=4, suppress=True):
            print(f" Proposed params are {x.T}")
        

        it += 1

        ## Graphing
        fig, ax = plt.subplots()
        ax.plot(OML_x,OML_z,color='black',label='Outer')
        booms, = plt.plot(table_calc[:,2],table_calc[:,3],\
            color='blue',linestyle=None,marker='o',label='Booms')
        ax.set_aspect('equal')
        booms.remove()
        booms, = plt.plot(table_calc[:,2],table_calc[:,3],\
                color='blue',linestyle=None,marker='o',label='Booms')
        plt.title('OML and booms')
        ax.set_aspect('equal')
        plt.legend()
        plt.show()

if __name__ == '__main__':
    iterate()
    with np.printoptions(precision=4, suppress=True):    
        print(f"*"*50)
        print(f"The final solution is {x.T}")
        print(f"The final minimum factors of safety are {G(x).T}")
        print(f"The final area is {wingbox_area(x)}m^2")
        print(f"*"*50)
        print(f"The initial area is {wingbox_area(init_vec)}m^2")
        print(f"The initial minimum factors of safety are {G(init_vec).T}")
        
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    
    ax1 = plt.subplot()
    for row in range(3):
        stress_MS, = ax1.plot(twenty_iter[row,:], color='blue')
    ax2 = ax1.twinx()
    for row in range(3,3+m):
        millimetres, = ax2.plot(twenty_iter[row,:], color='black')
    area, = ax2.plot(twenty_iter[-1,:], color='orange')
    
    ax1.set_ylabel('Stress margins of safety [dimensionless]')
    ax2.set_ylabel('Parameter dimensions [mm, mm2]')
    ax1.set_xlabel('Iteration')
    plt.legend([stress_MS, millimetres, area], ["MS", "Parameter dimensions", "Wingbox total area"])
    
    plt.xticks(np.arange(0,num_iter-1,2))
    plt.show()
    
    if write_to_excel == True:
        import pandas as pd
        
        filepath = 'final_stress_table_calculation.xlsx'
        
        column_names = ["i - Element no.",
                        "k - Element type",
                        "X - absolute x-coordinate",
                        "Z - absolute z-coordinate",
                        "x - x-coordinate relative to centroid",
                        "z - z-coordinate relative to centroid",
                        "Ae - Effective area",
                        "Ae*x^2",
                        "Ae*z^2",
                        "Ae*x*z",
                        "sigma_y - bending stresses",
                        "dx",
                        "dz",
                        "A_enclosed - triangular enclosed area of each panel section to centroid",
                        "q' - non-adjusted shear flow",
                        "q'*A_enclosed",
                        "q - true shear-flow",
                        "tau - shear stress",
                        "MS_tens",
                        "MS_crippling",
                        "MS_shear"
                        ]
        df = pd.DataFrame(table_calc)
        
        df.to_excel(filepath, index=False, header = column_names)
        
        print(f"*"*50)
        print(f"Stress calculation for final parameters has been output to file named {filepath}")
    
    pass
    
     


    



