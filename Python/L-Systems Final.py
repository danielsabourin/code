import math
import random
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import time

# PDD_item class holds the information of where the plotting "goes":
# i.e position (X,Y,Z)
# angles/facing (theta, phi)
# distance travelled = d (since it changes)
class PDD_item:
    def __init__(self,x,y,z,theta,phi,distance,width):
        self.X = x
        self.Y = y
        self.Z = z
        self.theta = theta
        self.phi = phi
        self.d = distance
        self.w = width 
    def __repr__(self): # for testing purposes
        return "(x,y,z): ({0.X},{0.Y},{0.Z}) \n theta: {0.theta} \n phi: {0.phi} \n distance: {0.d} \n width: {0.w} \n".format(self)

# the plot is global in order to be able to draw lines as function calls
fig = plt.figure()
ax = fig.gca(projection='3d')
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() 

# global constant of growth
golden_ratio = (1 + math.sqrt(5)) / 2

# copy_update(pdd, to_update, new_value) consumes a PDD_item
# and updates the field given by the string to_update, to the 
# value new_value
# copy_update: PDD_item, Str, (anyOf Float Int) -> PDD_item
def copy_update(pdd, to_update, new_value):
    new_pdd = PDD_item(pdd.X, pdd.Y, pdd.Z, pdd.theta, pdd.phi, pdd.d, pdd.w)
    if "X" in to_update: 
        new_pdd.X = new_value[new_value.index("X") + 1]
    if "Y" in to_update: 
        new_pdd.Y = new_value[new_value.index("Y") + 1]
    if "Z" in to_update: 
        new_pdd.Z = new_value[new_value.index("Z") + 1] 
    if "theta" in to_update: 
        new_pdd.theta = new_value[new_value.index("theta") + 1]   
    if "phi" in to_update: 
        new_pdd.phi = new_value[new_value.index("phi") + 1]   
    if "d" in to_update: 
        new_pdd.d = new_value[new_value.index("d") + 1]  
    if "w" in to_update: 
        new_pdd.w = new_value[new_value.index("w") + 1]      
    return new_pdd

# expand_L_system(axiom, rules, n) consumes an axiom as a String,
# rules as a dictionary of Strings, and an integer n:
# it produces one string which represents the L-system iterated n times using
# the axiom and rules found in the dictionary
# expand_L_system: String Dictionary Int -> String
def expand_L_system(axiom, rules, n):
    if n == 0:
        return axiom
    z = 0
    while z < len(axiom):
        if axiom[z] in rules:
            length = len(rules[axiom[z]])
            axiom = axiom[:z] + rules[axiom[z]] + axiom[z+1:]
            z = z + length
        else:
            z = z + 1
    return expand_L_system(axiom, rules, n-1)

# filter_placeholders(exp_L) consumes a the string exp_L and removes
# all of the placeholder characters
def filter_placeholders(exp_L):
    new_exp_L = ""
    allowed_list = ["F", "+", "-", "<", ">","[","]"]
    for letter in exp_L:
        if not (letter in allowed_list):
            pass
        else:
            new_exp_L = new_exp_L + letter
    return new_exp_L

# draws the actual branches
def draw_line(currPDD, d_f, multi_thresh, g_thresh):
    d_short = np.linspace(0, currPDD.d, 2)
    phi = currPDD.phi/360 * 2 * np.pi
    theta = currPDD.theta/360 * 2 * np.pi
    
    x_short = currPDD.X + (d_short * np.sin(phi) * np.cos(theta))
    y_short = currPDD.Y + (d_short * np.sin(phi) * np.sin(theta))
    z_short = currPDD.Z + (d_short * np.cos(phi))
    
    widths = []
    for i in range (0, 5):
        if i == 4:
            widths = widths + [currPDD.w - (i/5 * d_f) * 0.9]
        else:
            widths = widths + [currPDD.w - (i/5 * d_f) * (1 -(0.05 * i))]
    
    if currPDD.d < g_thresh:
        ax.plot(x_short,y_short,z_short, color = "green", linewidth = currPDD.w)
    elif currPDD.d > multi_thresh: # for larger widths, we use 5 small lengths to 
        # make the transition more aesthetic
        d_long = np.linspace(0, currPDD.d, 5)
        x_long = currPDD.X + (d_long * np.sin(phi) * np.cos(theta))
        y_long = currPDD.Y + (d_long * np.sin(phi) * np.sin(theta))
        z_long = currPDD.Z + (d_long * np.cos(phi)) 
        for i in range(0, 5):
            ax.plot(x_long[i:i+2],y_long[i:i+2],z_long[i:i+2],color="saddlebrown",linewidth=widths[i])
    else:
        ax.plot(x_short,y_short,z_short,color="saddlebrown", linewidth=currPDD.w)
    
# draw_L_system(axiom, rules, n, d, theta, phi) plots the L-system associated with
# the axiom, rules, number of iterations n, where
# d = distance travelled
# theta = angle (in degrees) to turn around the z-axis (think spherical co-ordinates)
# phi = angle (in degrees) to turn around the y-axis
# draw_L_system: String Dictionary Int Float Float Float -> None
# Effects:
# plots the L-system on the 3d grid
## NOTE ON ANGLES:
    ## + is counter clockwise around z-axis, - clockwise
    ## < is up along phi and > is down
    
def draw_L_system(exp_L, d, d_f, theta_L, phi_L, theta_i,pos_i,g_thresh): 
    
    # position (x,y,z), direction(theta, phi), theta_initial, distance(d), width(w) 
    PDD = [PDD_item(pos_i[0],pos_i[1],pos_i[2],theta_i, 0, d, d/2)] 
    bracket_index_list = []
    pdd_index = 0 # index of PDD list
    
    if d > 5:
        multi_thresh = d * 0.99 * (d_f ** 2)
    elif d > 2:
        multi_thresh = d * 0.99
    else:
        multi_thresh = d + 1    
    
    for letter in exp_L:
        currPDD = PDD[pdd_index]
        
        theta = random.normalvariate(theta_L[0],theta_L[1])
        phi = random.normalvariate(phi_L[0],phi_L[1])            
        
        if letter == "F":
            x = currPDD.X
            y = currPDD.Y
            z = currPDD.Z
            
            phi_c = currPDD.phi/360 * 2 * np.pi # phi in radians
            theta_c = currPDD.theta/360 * 2 * np.pi # theta in radians
            
            delta_x = currPDD.d * np.sin(phi_c) * np.cos(theta_c)
            x = x + delta_x
            delta_y = currPDD.d * np.sin(phi_c) * np.sin(theta_c)
            y = y + delta_y
            delta_z = currPDD.d * np.cos(phi_c)
            z = z + delta_z
            
            draw_line(currPDD, d_f, multi_thresh, g_thresh) # draws the line

            PDD = PDD + [copy_update(currPDD,["X", "Y", "Z"], ["X", x, "Y", y, "Z", z])]

        elif letter == "+":
            new_theta = (currPDD.theta + theta) % 360
            PDD = PDD + [copy_update(currPDD,["theta"], ["theta", new_theta])]
        
        elif letter == "-":
            if (currPDD.theta - theta) < 0:
                new_theta = (360 + (currPDD.theta - theta)) % 360 # keeps theta < 360
            else:
                new_theta = (currPDD.theta - theta) % 360 # keeps theta < 360
            PDD = PDD + [copy_update(currPDD,["theta"], ["theta", new_theta])]     
            
        elif letter == ">":
            if (currPDD.phi + phi) > 180: # when phi goes over 180, it rotates theta by 180
                new_phi = 180 - ((currPDD.phi + phi) % 180)
                new_theta = 360 - currPDD.theta
                PDD = PDD + [copy_update(currPDD, ["phi", "theta"], ["phi", new_phi, "theta", new_theta])]
            else:
                new_phi = currPDD.phi + phi
                PDD = PDD + [copy_update(currPDD, ["phi"], ["phi",new_phi])]   
                
        elif letter == "<":
            if (currPDD.phi - phi) < 0: # when phi goes under 0, it rotates theta by 180
                new_phi = abs(currPDD.phi - phi)
                new_theta = 360 - currPDD.theta
                PDD = PDD + [copy_update(currPDD, ["phi", "theta"], ["phi", new_phi, "theta", new_theta])]            
            else:
                new_phi = (currPDD.phi - phi)
                PDD = PDD + [copy_update(currPDD, ["phi"], ["phi",new_phi])] 
                
        elif letter == "[":
            bracket_index_list = bracket_index_list + [currPDD]
            PDD = PDD + [copy_update(currPDD, ["d", "w"], ["d", currPDD.d * d_f, "w", currPDD.w * d_f])]
            
        elif letter == "]":
            if len(bracket_index_list) == 1:
                PDD = PDD + [bracket_index_list[0]]
                bracket_index_list = []
            else:
                PDD = PDD + [bracket_index_list[-1]]
                bracket_index_list = bracket_index_list[:-1]
        pdd_index += 1
       

# draws a random forest based on the plants described by the list [trees]
# n is the number of iterations to run on the trees,
# x_space and y_space represent the spacing in between the trees and 
# x_dim, y_dim are the dimensions of the entire forest
def draw_forest(trees, n, x_space, y_space,x_dim, y_dim):
    x_max = int(math.floor((x_dim//x_space)))
    exp_L = []
    for elem in trees:
        exp_L = exp_L + [filter_placeholders(expand_L_system((elem[0])[0],(elem[0])[1],n))]
    
    for row in range(0, x_max):
        y_max = int(math.floor((y_dim//y_space)*random.uniform(1, 1.2)))
        for column in range(0, y_max):
            trees_random = random.uniform(-0.1, 99.9)
            tree_index = 0
            for elem in trees: # pick a random type of tree of the list given
                if trees_random < elem[1]:
                    L_sys = exp_L[tree_index]
                    d_L  = (elem[0])[2]
                    d_f  = (elem[0])[3]
                    theta  = (elem[0])[4]
                    phi  = (elem[0])[5]
                    g_thresh = (elem[0])[6]
                else:
                    trees_random = trees_random - elem[1]
                tree_index += 1
            
            d = random.normalvariate(d_L[0], d_L[1]) # randomized starting height
            d = max(d, d_L[0] - 3*d_L[1])  
                                         
            theta_i = random.normalvariate(180,60) # randomly rotates plants
            theta_i = max(0,theta_i)
            
            # trees in the forest are slightly randomly shifted to appear more natural
            tree_x = random.uniform(0, x_dim*0.2) + row * x_space + ((random.uniform(-1, 1)) * x_space * 0.3)
            tree_y = random.uniform(0, y_dim*0.2) + column * y_space + ((random.uniform(-1, 1)) * y_space * 0.3)   
            draw_L_system(L_sys, d, d_f, theta, phi, theta_i, [tree_x,tree_y,0],g_thresh)
      
    plt.show()   
    
# a plant is created as follows:
# plant = ["axiom", {productions}, distance(mean,std.dev), d_f, theta(mean,std.dev), phi(mean,std.dev), green_threshold)
# note on green threshold: if the length of a branch goes beneath this number, the branch will be green rather than brown
tree1 = ["FX", {"X": "[>---FX][>++++FX][>+++FX]"}, [10,2], 1/golden_ratio, [22,4], [12,2], 1]
broken_tree = ["FX", {"X": "[>+FX]"}, [9,1], 1/golden_ratio, [180,30], [3,0.5], -5]
shrub = ["X", {"X": "F[>--FX][>+++FX][>----FX][>+++++FX]"}, [4,1], 1/golden_ratio, [20,3], [10,2], 1]
grass = ["X", {"X": "[F>F>F>F>FX]+[F>F>FX]++[F>F>F>FX]+++[F>F>F>F>FX]--+++[F>F>FX]++++[F>F>FX]+++++[F>F>F>FX]"}, [1,0.2], 0.6, [35,5], [7,1], 5]

start_time = time.time()
draw_forest([[grass,100]], 2, 1.5, 1.5, 22, 22)
draw_forest([[shrub,40],[tree1,50],[broken_tree,10]], 6, 8, 8, 25, 25)
print("--- %s seconds ---" % (round(time.time() - start_time))) # just keeps track of how long the program takes to run