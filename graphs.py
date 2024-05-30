#%% IMPORTS
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

#%% DATA
# NOTE: length values (mm) | mass (grams) | diameter (mm)

df = pd.read_csv('dataset.csv') # creating a dataframe from the dataset.csv file
unique_fins = df['num_fins'].unique() # identifies the unique values in the 'num_fins' column
dataframes = {num_fins: df[df['num_fins'] == num_fins] for num_fins in unique_fins} # creates a dict. to hold dataframes

# creating dataframes for each set of different number of fins
df_500 = dataframes[500]
df_1k = dataframes[1000]
df_5k = dataframes[5000]
df_10k = dataframes[10000]

# sample length values for plotting
min_length = 8 # mm
max_length = 100 # mm
length_values = np.linspace(min_length, max_length, 200) # generating the values
length_values_5k = np.linspace(min_length-4, max_length, 200) # generating the values
length_values_10k = np.linspace(min_length-7.5, max_length, 200) # generating the values

# sample diameter values for plotting
min_diameter = 0.01 # mm
max_diameter = 10 # mm
diameter_values_500 = np.linspace(min_diameter, max_diameter, 200) # generating the values
diameter_values_1k = np.linspace(min_diameter, max_diameter-6, 200) # generating the values
diameter_values_5k = np.linspace(min_diameter, max_diameter-8, 200) # generating the values
diameter_values_10k = np.linspace(min_diameter-0.005, max_diameter-8, 200) # generating the values

#%% GRAPHS (500 FINS)
p500_L = np.polyfit(df_500['length'], df_500['mass'], 6) # creating a regression model of 500 fins length data
p500_D = np.polyfit(df_500['diameter'], df_500['mass'], 6) # creating a regression model of 500 fins diameter data

def fins_500_length(L): # creating a function/equation based on regression of length for plotting
    return np.polyval(p500_L, L)
def fins_500_diameter(D): # creating a function/equation based on regression of diameter for plotting
    return np.polyval(p500_D, D)

plt.plot(length_values, fins_500_length(length_values), label='Length') # generating lineplot for length data
plt.plot(diameter_values_500, fins_500_diameter(diameter_values_500), label='Diameter') # generating lineplot for mass data
plt.scatter(df_500['length'], df_500['mass'], marker='x') # generating scatter plot from real data points on length
plt.scatter(df_500['diameter'], df_500['mass'], marker='o') # generating scatter plot from real data points on diameter

# decorating the plot
plt.grid()
plt.legend(title='Dimension of Fin' )
plt.xlabel('Length & Diameter (mm)')
plt.ylabel('Mass (g)')
plt.title('Length of Fin & Fin Diameter vs. Total Mass of Fins for 500 Fins')
plt.show()

#%% GRAPHS (1000 FINS)
p1k_L = np.polyfit(df_1k['length'], df_1k['mass'], 6) # creating a regression model of 500 fins length data
p1k_D = np.polyfit(df_1k['diameter'], df_1k['mass'], 6) # creating a regression model of 500 fins diameter data

def fins_1k_length(L): # creating a function/equation based on regression of length for plotting
    return np.polyval(p1k_L, L)
def fins_1k_diameter(D): # creating a function/equation based on regression of diameter for plotting
    return np.polyval(p1k_D, D)

plt.plot(length_values, fins_1k_length(length_values), label='Length') # generating lineplot for length data
plt.plot(diameter_values_1k, fins_1k_diameter(diameter_values_1k), label='Diameter') # generating lineplot for mass data
plt.scatter(df_1k['length'], df_1k['mass'], marker='x') # generating scatter plot from real data points on length
plt.scatter(df_1k['diameter'], df_1k['mass'], marker='o') # generating scatter plot from real data points on diameter

# decorating the plot
plt.grid()
plt.legend(title='Dimension of Fin' )
plt.xlabel('Length & Diameter (mm)')
plt.ylabel('Mass (g)')
plt.title('Length of Fin & Fin Diameter vs. Total Mass of Fins for 1000 Fins')
plt.show()

#%% GRAPHS (5000 FINS)
p5k_L = np.polyfit(df_5k['length'], df_5k['mass'], 6) # creating a regression model of 500 fins length data
p5k_D = np.polyfit(df_5k['diameter'], df_5k['mass'], 6) # creating a regression model of 500 fins diameter data

def fins_5k_length(L): # creating a function/equation based on regression of length for plotting
    return np.polyval(p5k_L, L)
def fins_5k_diameter(D): # creating a function/equation based on regression of diameter for plotting
    return np.polyval(p5k_D, D)

plt.plot(length_values_5k, fins_5k_length(length_values_5k), label='Length') # generating lineplot for length data
plt.plot(diameter_values_5k, fins_5k_diameter(diameter_values_5k), label='Diameter') # generating lineplot for mass data
plt.scatter(df_5k['length'], df_5k['mass'], marker='x') # generating scatter plot from real data points on length
plt.scatter(df_5k['diameter'], df_5k['mass'], marker='o') # generating scatter plot from real data points on diameter

# decorating the plot
plt.grid()
plt.legend(title='Dimension of Fin' )
plt.xlabel('Length & Diameter (mm)')
plt.ylabel('Mass (g)')
plt.title('Length of Fin & Fin Diameter vs. Total Mass of Fins for 5000 Fins')
plt.show()

#%% GRAPHS (10000 FINS)
p10k_L = np.polyfit(df_10k['length'], df_10k['mass'], 2) # creating a regression model of 500 fins length data
p10k_D = np.polyfit(df_10k['diameter'], df_10k['mass'], 2) # creating a regression model of 500 fins diameter data

def fins_10k_length(L): # creating a function/equation based on regression of length for plotting
    return np.polyval(p10k_L, L)
def fins_10k_diameter(D): # creating a function/equation based on regression of diameter for plotting
    return np.polyval(p10k_D, D)

plt.plot(length_values_10k, fins_10k_length(length_values_10k), label='Length') # generating lineplot for length data
plt.plot(diameter_values_10k, fins_10k_diameter(diameter_values_10k), label='Diameter') # generating lineplot for mass data
plt.scatter(df_10k['length'], df_10k['mass'], marker='x') # generating scatter plot from real data points on length
plt.scatter(df_10k['diameter'], df_10k['mass'], marker='o') # generating scatter plot from real data points on diameter

# decorating the plot
plt.grid()
plt.legend(title='Dimension of Fin' )
plt.xlabel('Length & Diameter (mm)')
plt.ylabel('Mass (g)')
plt.title('Length of Fin & Fin Diameter vs. Total Mass of Fins for 10000 Fins')
plt.show()

#%% ALL GRAPHS TOGETHER
plt.plot(length_values, fins_500_length(length_values), label='Length (500 fins)') # generating lineplot for length data
plt.plot(diameter_values_500, fins_500_diameter(diameter_values_500), label='Diameter (500 fins)') # generating lineplot for mass data

plt.plot(length_values, fins_1k_length(length_values), label='Length (1k fins)') # generating lineplot for length data
plt.plot(diameter_values_1k, fins_1k_diameter(diameter_values_1k), label='Diameter (1k fins)') # generating lineplot for mass data

plt.plot(length_values_5k, fins_5k_length(length_values_5k), label='Length (5k fins)') # generating lineplot for length data
plt.plot(diameter_values_5k, fins_5k_diameter(diameter_values_5k), label='Diameter (5k fins)') # generating lineplot for mass data

plt.plot(length_values_10k, fins_10k_length(length_values_10k), label='Length (10k fins)') # generating lineplot for length data
plt.plot(diameter_values_10k, fins_10k_diameter(diameter_values_10k), label='Diameter (10k fins)') # generating lineplot for mass data

plt.title('Plot of All Fin Numbers of Length & Diameter vs. Mass')
plt.xlabel('Fin Diameter & Length (mm)')
plt.ylabel('Total Fin Mass (g)')
plt.grid()
plt.legend()
plt.show()