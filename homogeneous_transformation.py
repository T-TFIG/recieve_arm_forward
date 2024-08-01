#!/usr/bin/env python3
import numpy as np
import sympy as sp

def dh_transform(theta, d, a, alpha):
    """
    Computes the homogeneous transformation matrix using DH parameters.
    
    Parameters:
    theta (float): Joint angle (radians)
    d (float): Link offset
    a (float): Link length
    alpha (float): Link twist (radians)
    
    Returns:
    sympy.Matrix: 4x4 homogeneous transformation matrix
    """
    # Create the transformation matrix
    T = sp.Matrix([
        [sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
        [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
        [0,              sp.sin(alpha),               sp.cos(alpha),               d],
        [0,              0,                           0,                           1]
    ])
    
    return sp.simplify(T)

def transform_frame_to_end_effector(dh_table):
    """
    Computes the transformation matrix from the base frame to the end effector using a DH table.
    
    Parameters:
    dh_table (numpy.ndarray): DH table where each row represents [theta, d, a, alpha]
    
    Returns:
    sympy.Matrix: The overall transformation matrix
    """
    T = sp.eye(4)  # Initialize as identity matrix
    for i in range(dh_table.shape[0]):
        theta, d, a, alpha = dh_table[i]
        T_i = dh_transform(theta, d, a, alpha)
        T = sp.simplify(T * T_i)  # Simplify at each step to reduce complexity
    
    return T

# Define symbolic variables for DH parameters
l1, l2, l3, l4, l5, l6, l7 = sp.symbols('l1 l2 l3 l4 l5 l6 l7')
zeta_1, zeta_2, zeta_3, zeta_4 = sp.symbols('zeta_1 zeta_2 zeta_3 zeta_4')

# Define the DH table as a numpy array
dh_table = np.array([
    [0, sp.pi/2, 0, sp.pi/2],
    [l1, 0, 0, zeta_1],
    [0, -sp.pi/2, -l2, -sp.pi/2],
    [l3, 0, 0, zeta_2 + sp.pi/2],
    [l4, 0, l5, -sp.pi/2],
    [0, sp.pi/2, 0, zeta_3],
    [l6, 0, 0, zeta_4],
    [l7, 0, 0, 0]
])

# Compute the transformation matrix from base frame to end effector
T_end_effector = transform_frame_to_end_effector(dh_table)

# Simplify the x, y, z positions
x_position = sp.simplify(T_end_effector[0, 3])
y_position = sp.simplify(T_end_effector[1, 3])
z_position = sp.simplify(T_end_effector[2, 3])

print("X Position:")
print(x_position)
print("Y Position:")
print(y_position)
print("Z Position:")
print(z_position)




# X Position:
# -l2*cos(l1) - l5*(sin(l1)*sin(l3)*sin(zeta_1) - cos(l1)*cos(l3))*cos(l4) + l5*(sin(l1)*sin(zeta_1)*sin(zeta_2)*cos(l3) 
# - sin(l1)*cos(zeta_1)*cos(zeta_2) + sin(l3)*sin(zeta_2)*cos(l1))*sin(l4) + pi*((sin(l1)*sin(l3)*sin(zeta_1) - cos(l1)*cos(l3))*sin(l4) 
# + (sin(l1)*sin(zeta_1)*sin(zeta_2)*cos(l3) 
# - sin(l1)*cos(zeta_1)*cos(zeta_2) + sin(l3)*sin(zeta_2)*cos(l1))*cos(l4))/2 - pi*sin(l1)*sin(zeta_1)/2


# Y Position:
# -l5*(sin(zeta_1)*cos(zeta_2) + sin(zeta_2)*cos(l3)*cos(zeta_1))*sin(l4) + l5*sin(l3)*cos(l4)*cos(zeta_1) 
# - pi*((sin(zeta_1)*cos(zeta_2) + sin(zeta_2)*cos(l3)*cos(zeta_1))*cos(l4) + sin(l3)*sin(l4)*cos(zeta_1))/2 + pi*cos(zeta_1)/2


# Z Position:
# -l2*sin(l1) + l5*(sin(l1)*cos(l3) + sin(l3)*sin(zeta_1)*cos(l1))*cos(l4) + l5*(sin(l1)*sin(l3)*sin(zeta_2) 
# - sin(zeta_1)*sin(zeta_2)*cos(l1)*cos(l3) + cos(l1)*cos(zeta_1)*cos(zeta_2))*sin(l4) 
# - pi*((sin(l1)*cos(l3) + sin(l3)*sin(zeta_1)*cos(l1))*sin(l4) - (sin(l1)*sin(l3)*sin(zeta_2) 
# - sin(zeta_1)*sin(zeta_2)*cos(l1)*cos(l3) + cos(l1)*cos(zeta_1)*cos(zeta_2))*cos(l4))/2 + pi*sin(zeta_1)*cos(l1)/2 + pi/2