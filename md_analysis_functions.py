"""
Module contains functions used to analyse data from CASTEP .md files
"""
import numpy as np
import scipy.stats as stats
import math
import statistics


# contains functions used to analyse data from castep .md files
def get_data(atom_list, md_file, value):
    """
    Reads the .md file and returns a specified set of vector data
    Applies to every atom specified in the atom_list
    Also returns an atom identifier 'ID' based on the atomic symbol and the number

    Parameters
    ----------
    atom_list : list
        Atoms of interest to get the data for
        e.g. ["H1", "H2", "H3"]
        Can select atoms using the atomic symbol and atom number as appears in the .md file
        e.g. "H                 1"  has an ID "H1"

    md_file : string
        The file to read

    value : string
        The data set to return
        'R' Position, 'V' Velocity, 'F' Force
        If letters are present due to atomic symbols can use '<-- R' Position, '<-- V' Velocity, '<-- F' Force

    Returns
    -------
    data : list of arrays
        Data for each atom at every time step
        In the format data[atom][time][Atom, AtomNumber, x,y,z, value, ID]
    """

    # create empty container for data
    data = []
    # loop through all the provided atoms
    for atom in atom_list:
        # read the provided .md file
        with open(md_file, "r") as f:
            # create empty container for every timestep data
            line_lists = []
            # evaluate every line of the file
            for line in f:
                # if the quantity is the one asked for
                if value in line:
                    # split up the line removing spaces
                    line_list = line.split()
                    # create identifier based on the atomic symbol and atom number
                    identifier = line_list[0] + line_list[1]
                    # if the identifier matches the provided atom string
                    if atom == identifier:
                        # add the identifier to the line data
                        line_list.append(identifier)
                        # append the data for this timestep to the list of all timestamps
                        line_lists.append(line_list)
            # append the data for each atom to a single list
            data.append(line_lists)
    return data


def additional_points(atoms, time, md_file, new_point_dis, new_point_name):
    """
    Reads the .md file generates new points in the x y plain relative to specified atoms for a specified timestep

    Parameters
    ----------
    atoms : list
        Atoms of interest selected by ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    time : integer
        time step of interest

    md_file : string
        The file to read

    new_point_dis : list of lists
        distance of new points from existing point in Angstroms
        each item in the list consists of an x and y distance [[x1, y1], [x2, y2]]

    new_point_name : list
        corresponding identifiers for the new points
        e.g. ["position1", "position2"]

    Returns
    -------
    new_x_point_list : list
        new x coordinates

    new_y_point_list : list
        new y coordinates

    new_point_name_list : list
        identifiers for each new point
    """
    # create empty containers for coordinates and identifiers
    new_x_point_list = []
    new_y_point_list = []
    new_point_name_list = []

    # get the positional data for the selected atoms and loop
    data = get_data(atoms, md_file, "R")

    for i in range(0, len(data), 1):
        for j in range(0, len(new_point_dis)):
            # positional data from md files are in a.u
            # convert to angstroms here by multiplying by 0.529177249
            # then calculate the new positions
            new_x = float(data[i][time][2]) * 0.529177249 + new_point_dis[j][0]
            new_y = float(data[i][time][3]) * 0.529177249 + new_point_dis[j][1]
            # append the new coordinates to lists of x and y coordinates
            new_x_point_list.append(new_x)
            new_y_point_list.append(new_y)
            # generate a new coordinate identifier
            new_point_name_list.append(new_point_name[j])

    return new_x_point_list, new_y_point_list, new_point_name_list


def atom_distances(md_file, relations):
    """
     Calculate the distance between specified pairs of atoms


     Parameters
     ----------
     relations : list of lists
         Specify the atoms to find the distance between
         Each item a list of 2 values e.g. [['H1', 'H2'], ['H1', 'O1'], ['H2', 'O1']]

     md_file : string
         The file to read

     Returns
     -------
     distances_list : list of lists
         Distances between each relation for every time step in Angstroms
         distances_list[relation][time][distance]

     relation_list : list
         Order of corresponding relation IDs e.g. [['H1-H2'], ['H2-H3']]

     x_y_z_distances : list of lists
     """

    # container to hold the distances at each time step for each pair of atoms
    distances_list = []
    # relationship identifiers that correspond to the distances list
    relation_list = []
    # container to hold the separate x y and z distances
    x_y_z_distances_list = []

    # container to hold the atom identifiers in the relation
    atoms = []
    for relation in relations:
        atoms.append(relation[0])
        atoms.append(relation[1])
    # list of atoms to get the positional data for
    atom_list = list(set(atoms))
    data = get_data(atom_list, md_file, 'R')

    # For every specified pair
    for relation in relations:
        # create a unique identifier for the pair and store
        relation_id = relation[0] + '-' + relation[1]
        relation_list.append(relation_id)

        # Loop though the data to find the data associated with each atom in the pair
        # Remember data is in the format data[atom][time][Atom, AtomNumber, x,y,z, value, ID]
        for atom1 in range(0, len(data), 1):
            # if the string of the first atom matches a line of data
            # the integer values of atom1 gives the position of the specified atom data in the data list
            if relation[0] == data[atom1][0][7]:
                # do the same to find the corresponding integer value for the second atom
                for atom2 in range(0, len(data), 1):
                    if relation[1] == data[atom2][0][7]:

                        # with the atom1 and atom2 integers for data[atom1] and data[atom2]
                        # these can be used to calculate the distance between their coordinates
                        # containers of distances for the current pair of atoms at every timestep
                        distances = []
                        x_y_z_distances = []
                        # Loop through every time step
                        for time in range(0, len(data[0]), 1):
                            # subtract x y and z coordinates
                            x_distance = float(data[atom1][time][2]) - float(data[atom2][time][2])
                            y_distance = float(data[atom1][time][3]) - float(data[atom2][time][3])
                            z_distance = float(data[atom1][time][4]) - float(data[atom2][time][4])
                            # append the separate x y and z coordinates to a master list
                            x_y_z_distances.append([x_distance * 0.529177249,
                                                    y_distance * 0.529177249,
                                                    z_distance * 0.529177249])
                            # find the distance in 3D using pythagoras
                            distance = math.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
                            # positions are in a.u so converting to Angstroms here
                            distances.append(distance * 0.529177249)

        # having calculated the distances of atoms for a particular relation append this to a combined lists
        distances_list.append(distances)
        x_y_z_distances_list.append(x_y_z_distances)

    # distances_list[relation][time][distance] relation_list[relation]
    return distances_list, relation_list, x_y_z_distances_list


def center_of_mass(atoms, md_file, masses):
    """
    gives the coordinates of the center of mass for a molecule given atom positions

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as generated from the 'get_data' function e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    x_co : list
        center of mass x coordinate for every time step in Angstroms

    y_co : list
        center of mass y coordinate for every time step in Angstroms

    z_co : list
        center of mass z coordinate for every time step in Angstroms
    """

    # get the position data of the atoms
    data = get_data(atoms, md_file, 'R')
    # calculate total mass of molecule
    total_mass = 0
    for mass in masses:
        total_mass = total_mass + mass

    # containers for the center of mass coordinate along each axis
    x_co = []
    y_co = []
    z_co = []
    # loop through each timestep of data
    for t in range(0, len(data[0]), 1):
        center_x = 0
        center_y = 0
        center_z = 0
        # summation of position * mass followed division of total mass for all atoms to give center of mass
        for i in range(0, len(atoms), 1):
            # data[atom][time][[Atom, AtomNumber, x,y,z, value, ID]]
            center_x = center_x + (float(data[i][t][2]) * masses[i])
            center_y = center_y + (float(data[i][t][3]) * masses[i])
            center_z = center_z + (float(data[i][t][4]) * masses[i])
        # convert a.u to angstroms
        x_co.append((center_x * 0.529177249) / total_mass)
        y_co.append((center_y * 0.529177249) / total_mass)
        z_co.append((center_z * 0.529177249) / total_mass)

    # return lists containing center of mass for each time step along each axis
    return x_co, y_co, z_co


def center_of_velocity(atoms, md_file, masses):
    """
    gives the coordinates of the center of velocity for a molecule

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    x_co : list
        center of mass x coordinate for every time step in a.u

    y_co : list
        center of mass y coordinate for every time step in a.u

    z_co : list
        center of mass z coordinate for every time step in a.u
    """

    # get the velocity data of the atoms
    data = get_data(atoms, md_file, 'V')
    # calculate total mass of molecule
    total_mass = 0
    for mass in masses:
        total_mass = total_mass + mass

    # containers for the center of velocity along each axis
    x_co = []
    y_co = []
    z_co = []
    # loop through each timestep of data
    for t in range(0, len(data[0]), 1):
        center_x = 0
        center_y = 0
        center_z = 0
        # summation of velocity * mass followed division of total mass for all atoms to give center of mass
        for i in range(0, len(atoms), 1):
            # data[atom][time][[Atom, AtomNumber, x,y,z, value, ID]]
            center_x = center_x + (float(data[i][t][2]) * masses[i])
            center_y = center_y + (float(data[i][t][3]) * masses[i])
            center_z = center_z + (float(data[i][t][4]) * masses[i])
        x_co.append(center_x / total_mass)
        y_co.append(center_y / total_mass)
        z_co.append(center_z / total_mass)

    # return lists containing center of velocity for each time step along each axis
    return x_co, y_co, z_co


def find_closest(input_list, input_value):
    """
    finds the closest value in a list to a given input

    Parameters
    ----------
    input_list: list of floats
        list to search through

    input_value: float
        value to find closest of

    Returns
    -------
    output_value: float
        the closest value to the input_value that is in input_list
    """
    # set the input list as a list of floats
    arr = np.asarray(input_list)

    # subtract the input_value for every element in the array and find the absolute value
    difference_array = np.abs(arr - input_value)
    # find the index with the of the smallest value
    index = difference_array.argmin()
    # return the value of the index
    output_value = arr[index]

    return output_value


def get_ads_energy(md_file, mol_energy, sub_energy):
    """
    calculates the adsorption energy given the energy of a system and its separate components

    Parameters
    ----------
    md_file : string
         The file to read

    mol_energy : float
        Total energy of molecule without substrate (eV)

    sub_energy : float
        Total energy of substrate without molecule (eV)

    Returns
    -------
    adsorption_energy_array : list
        The adsorption energy for every time step (eV)
    """
    # calculate the energy of the molecule and substrate together
    reference_energy = mol_energy + sub_energy

    # get the energies at every timestep from the md file in eV
    time, total_energy, hamiltonian_energy, kinetic_energy = get_energy(md_file)

    # container to hold the adsorption energy of each timestep
    adsorption_energy_array = []

    # for every time step energy
    for i in range(0, len(total_energy)):
        # calculate adsorption energy (E_ads = E_ab - (E_a + E_b))
        e_ads = total_energy[i] - reference_energy
        # append the values to a list
        adsorption_energy_array.append(e_ads)

    return adsorption_energy_array


def angle_between(a, b):
    """
    Calculates angle between two vectors

    Parameters
    ----------
    a : list
        vector a

    b : list
        vector b

    Returns
    -------
    angle:
        angle between each vector in radians

    """
    # calculate dot product to give a*b*cos(angle)
    dot_prod = np.dot(a, b)
    # calculate cross product to give a*b*sine(angle)
    cross_prod = np.cross(a, b)
    # combine both and compute angle (cancel out magnitudes of a and b as tan(angle)=sin(angle)/cos(angle))
    angle = np.arctan2(cross_prod, dot_prod)

    return angle


def get_energy(md_file):
    """
    Retrieves the energy data from the .md file

    Parameters
    ----------
    md_file : string
         The file to read

    Returns
    -------
    time : list
        the number of time steps as a list

    total_energy : list
        the total energy for every time step

    hamiltonian_energy : list
        the hamiltonian_energy for every time step

    kinetic_energy : list
        the kinetic_energy for every time step
    """

    # container to hold the data
    data = []
    # read the md file

    with open(md_file, "r") as f:
        # loop round every line of the file
        for line in f:
            # .md files contain <-- E for the energy
            if '<-- E' in line:
                # split the data into a list removing white spaces
                line = line.split()
                data.append(line)

    # returns the number of timestamps
    time = list(range(0, len(data)))

    # for conversion of Hartree's to eV
    to_ev = 27.211324570273

    # containers for the different energies
    total_energy = []
    hamiltonian_energy = []
    kinetic_energy = []
    # loop round every timestep
    for i in range(0, len(data), 1):
        # get the energy for each line stored as data[timestep][total, hamiltonian, kinetic]
        # convert to eV
        total_energy.append(float(data[i][0]) * to_ev)
        hamiltonian_energy.append(float(data[i][1]) * to_ev)
        kinetic_energy.append(float(data[i][2]) * to_ev)

    # Time, Total, Hamiltonian, Kinetic energies
    return time, total_energy, hamiltonian_energy, kinetic_energy


def get_unit_cell(md_file):
    """
    Gets the unit cell parameters of the .md file
    (provided remain fixed for the simulation)

    Parameters
    ----------
    md_file : string
         The file to read

    Returns
    -------
    constants: list
        unit cell lattice constants a, b and c in Angstroms

     angles: list
        unit cell lattice constants alpha, beta, gamma in degrees
    """

    # container to hold lattice constant data
    data = []
    # read the file
    with open(md_file, "r") as f:
        # loop through every line
        for line in f:
            # lines with '<-- h' contain the unit cell data
            if '<-- h' in line:
                line = line.split()
                data.append(line)
            # only need the 3 lines of cartesian components (provided cell is fixed)
            if len(data) == 3:
                break

        # data should be a list of size 3 3 (plus the < -- h)

    # to convert a.u to angstroms
    to_ang = 0.529177210903

    # organise the cartesian components of lattice vectors into a single venerable and convert to angstroms
    a_cartesian = np.array([float(data[0][0]), float(data[0][1]), float(data[0][2])]) * to_ang
    b_cartesian = np.array([float(data[1][0]), float(data[1][1]), float(data[1][2])]) * to_ang
    c_cartesian = np.array([float(data[2][0]), float(data[2][1]), float(data[2][2])]) * to_ang

    # Calculate the magnitude of the cartesian components to give the Lattice constants (Pythagorean)
    a = np.linalg.norm(a_cartesian)
    b = np.linalg.norm(b_cartesian)
    c = np.linalg.norm(c_cartesian)

    # Get the cos of each lattice angle by finding dot products of the cartesian / magnitude of two vectors
    cos_alpha = np.dot(b_cartesian, c_cartesian) / (b * c)  # alpha is the angle between b and c etc...
    cos_beta = np.dot(a_cartesian, c_cartesian) / (a * c)
    cos_gamma = np.dot(a_cartesian, b_cartesian) / (a * b)

    # use arc_cos to find each of the angles
    alpha = np.arccos(cos_alpha)
    beta = np.arccos(cos_beta)
    gamma = np.arccos(cos_gamma)

    # Convert angles to degrees
    alpha_deg = np.degrees(alpha)
    beta_deg = np.degrees(beta)
    gamma_deg = np.degrees(gamma)

    # arrange in lists
    constants = [a, b, c]
    angles = [alpha_deg, beta_deg, gamma_deg]

    return constants, angles


def get_furthest_distance(atoms, md_file, length):
    """
    Get total distance traveled in simulation (xy plane)

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    length : integer
        Total number of time steps in the simulation

    Returns
    -------
    xy_distance_time : float
        distance traveled from the starting to ending time step

    sum_of_distances : float
        the total distance traveled in all time steps
    """

    # get position data of atoms
    data = get_data(atoms, md_file, "R")

    # find total distance traveled and convert distances to Angstroms
    x_distance = abs(float(data[0][0][2]) - float(data[0][length][2])) * 0.529177249
    y_distance = abs(float(data[0][0][3]) - float(data[0][length][3])) * 0.529177249
    xy_distance = math.sqrt(x_distance ** 2 + y_distance ** 2)
    # distance traveled over course of the simulation
    xy_distance_time = xy_distance

    # caculate distance traveled in every timestep to give total distance traveled
    sum_of_distances = 0
    for t in range(0, length - 1):
        x_distance1 = abs(float(data[0][t][2]) - float(data[0][t + 1][2])) * 0.529177249
        y_distance1 = abs(float(data[0][t][3]) - float(data[0][t + 1][3])) * 0.529177249
        xy_distance1 = math.sqrt(x_distance1 ** 2 + y_distance1 ** 2)
        sum_of_distances = sum_of_distances + xy_distance1

    return xy_distance_time, sum_of_distances


def calculate_distances_from_start(atoms, md_file, masses):
    """
    Calculate the distance from the start for each timestep and the total distance traveled.

    Parameters:
    atoms : int
        Number of atoms in the molecule
    md_file : str
        Path to the file containing molecular dynamics data

    Returns:
    - distances_from_start: Array of distances from the start for each timestep
    - sum_of_distances: Total distance traveled over all timesteps
    """

    # get center of mass
    x_cos, y_cos, z_cos = center_of_mass(atoms, md_file, masses)

    # Extract x and y coordinates
    x_coords = np.asarray(x_cos)
    y_coords = np.asarray(y_cos)
    z_coords = np.asarray(z_cos)

    # Start position at time step 0
    start_position = np.array([x_coords[0], y_coords[0], z_coords[0]]).T

    # Calculate distances from the start for each timestep
    distances_from_start = np.sqrt((x_coords - start_position[0]) ** 2
                                   + (y_coords - start_position[1]) ** 2
                                   + (z_coords - start_position[2]) ** 2)

    # Extract distances for each timestep
    distances_from_start_per_timestep = distances_from_start

    # Calculate distances traveled between timesteps
    delta_x = np.diff(x_coords, axis=0)
    delta_y = np.diff(y_coords, axis=0)
    delta_z = np.diff(z_coords, axis=0)
    distances_traveled = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)

    # Sum the distances traveled
    sum_of_distances = np.sum(distances_traveled)
    return distances_from_start_per_timestep, sum_of_distances


def get_step_direction(atoms, md_file, masses):
    """
    Get a vector for the direction of motion for each time step

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    x_vector_list : list
        value of x vector for every time step

    y_vector_list : list
        value of y vector for every time step

    y_vector_list : list
        value of z vector for every time step
    """
    # find the center of mass for the specified atoms at each time step in angstroms
    x_cos, y_cos, z_cos = center_of_mass(atoms, md_file, masses)

    # containers for holding the x, y, z direction vectors for each time step
    x_vector_list = []
    y_vector_list = []
    z_vector_list = []
    # loop through all the timesteps (except the last as no new direction data is added)
    for i in range(0, len(x_cos) - 1):
        # the direction can be expressed by the magnitude of the difference in coordinate
        x_vector = x_cos[i + 1] - x_cos[i]
        # append the difference
        x_vector_list.append(x_vector)
        # repeat for y and z
        y_vector = y_cos[i + 1] - y_cos[i]
        y_vector_list.append(y_vector)
        z_vector = z_cos[i + 1] - z_cos[i]
        z_vector_list.append(z_vector)

    return x_vector_list, y_vector_list, z_vector_list


def get_temperature(md_file):
    """
    Get temperatures for each time step

    Parameters
    ----------
    md_file : string
        The file to read

    Returns
    -------
    temperatures : list
        temperature fore each time step

    """
    # container to hold temperature line data
    data = []
    # read the file
    with open(md_file, "r") as f:
        # loop through each line
        for line in f:
            # temperature line contains <-- T
            if '<-- T' in line:
                # remove spaces and append
                line = line.split()
                data.append(line)

    # container for the temperature at every time step
    temperatures = []

    # for every line that contains a temperature
    for i in range(0, len(data), 1):
        # the temperature is at index 0
        # convert from a.u to K
        temperatures.append(float(data[i][0]) * 3.1577464E5)

    return temperatures


def get_velocities(atoms, md_file, masses):
    """
    Get velocities for atoms and molecules at each time step

    Parameters
    ----------

    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    atom_xy_vs_list : list of arrays
        the velocity in the xy plain for every atom for every time step a.u
        atom_xy_vs_list[atom][time]

    atom_xyz_vs_list : list of arrays
        the velocity in the 3D for every atom for every time step a.u
        atom_xyz_vs_list[atom][time]

    mol_xy_list : list
        the velocity in the xy plain for the molecule every time step a.u

    mol_xyx_array : list
        the velocity in the xyz plain for the molecule every time step a.u
    """
    # get velocity data for selected atoms from file
    # data[atom][time][Atom, AtomNumber, x,y,z]
    data = get_data(atoms, md_file, 'V')

    # container for xy plain velocity of each atom
    atom_xy_vs_list = []
    # container for xyz plain velocity of each atom
    atom_xyz_vs_list = []
    # container for xy plain velocity of the molecule
    mol_xy_list = []
    # container for xy plain velocity of the molecule
    mol_xyz_list = []

    # loop over each atom
    for a in range(0, len(data)):
        # containers for atom velocities
        atom_xy_vs = []
        atom_xyz_vs = []
        # loop over evey timestep
        for t in range(0, len(data[0])):
            # get velocity data data[atom][time][Atom, AtomNumber, x,y,z]
            x_velocity = float(data[a][t][2])
            y_velocity = float(data[a][t][3])
            z_velocity = float(data[a][t][4])

            # calculate magnitudes in each plain
            xy_velocity = (math.sqrt(x_velocity ** 2 + y_velocity ** 2))
            xyz_velocity = (math.sqrt(x_velocity ** 2 + y_velocity ** 2 + z_velocity ** 2))

            # append the atom velocities to a list
            atom_xy_vs.append(xy_velocity)
            atom_xyz_vs.append(xyz_velocity)

        # append the velocities of every time step for the atom (atom_xy_vs_list[atom][time])
        atom_xy_vs_list.append(atom_xy_vs)
        atom_xyz_vs_list.append(atom_xyz_vs)

    total_mass = sum(masses)
    # for every time step
    for t in range(0, len(atom_xyz_vs_list[0])):
        # initialise velocities
        mol_xy = 0
        mol_xyz = 0
        # for every atom
        for a in range(0, len(atom_xyz_vs_list)):
            # find the velocity of the molecule as a whole
            # (sum of (position * masses)) / total mass
            mol_xy = mol_xy + atom_xy_vs_list[a][t] * masses[a]
            mol_xyz = mol_xyz + atom_xyz_vs_list[a][t] * masses[a]

        mol_xy_list.append(mol_xy / total_mass)
        mol_xyz_list.append(mol_xyz / total_mass)

    # atom list [atom][time]
    # molecule list [time]
    return atom_xy_vs_list, atom_xyz_vs_list, mol_xy_list, mol_xyz_list


def get_forces(atoms, md_file, masses):
    """
    Get forces for atoms and molecules at each time step

    Parameters
    ----------

    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    atom_xy_vs_list : list of arrays
        the force in the xy plain for every atom for every time step a.u
        atom_xy_vs_list[atom][time]

    atom_xyz_vs_list : list of arrays
        the force in the 3D for every atom for every time step a.u
        atom_xyz_vs_list[atom][time]

    mol_xy_list : list
        the force in the xy plain for the molecule every time step a.u

    mol_xyx_array : list
        the force in the xyz plain for the molecule every time step a.u
    """

    # get force data from file
    data = get_data(atoms, md_file, 'F')

    # containers for forces for each atom at every timestep
    atom_xy_vs_list = []
    atom_xyz_vs_list = []
    mol_xy_list = []
    mol_xyz_list = []
    # loop through each atom
    for a in range(0, len(data)):
        # containers for forces at every timestep
        atom_xy_vs = []
        atom_xyz_vs = []
        # loop round every timestep
        for t in range(0, len(data[0])):
            # get force data data[atom][time][Atom, AtomNumber, x,y,z]
            x_force = float(data[a][t][2])
            y_force = float(data[a][t][3])
            z_force = float(data[a][t][4])

            # calculate magnitudes in each plain
            xy_force = (math.sqrt(x_force ** 2 + y_force ** 2))
            xyz_force = (math.sqrt(x_force ** 2 + y_force ** 2 + z_force ** 2))

            # append the atom forces to a list
            atom_xy_vs.append(xy_force)
            atom_xyz_vs.append(xyz_force)

        # append the velocities of every time step for the atom (atom_xy_vs_list[atom][time])
        atom_xy_vs_list.append(atom_xy_vs)
        atom_xyz_vs_list.append(atom_xyz_vs)

    total_mass = sum(masses)
    # for every time step
    for t in range(0, len(atom_xyz_vs_list[0])):
        # initialise forces
        mol_xy = 0
        mol_xyz = 0
        # for every atom
        for a in range(0, len(atom_xyz_vs_list)):
            # find the force of the molecule as a whole
            # (sum of (position * masses)) / total mass
            mol_xy = mol_xy + atom_xy_vs_list[a][t] * masses[a]
            mol_xyz = mol_xyz + atom_xyz_vs_list[a][t] * masses[a]

        mol_xy_list.append(mol_xy / total_mass)
        mol_xyz_list.append(mol_xyz / total_mass)

    # atom list [atom][time]
    # molecule list [time]
    return atom_xy_vs_list, atom_xyz_vs_list, mol_xy_list, mol_xyz_list


def get_translational_energy(velocities, masses):
    """
    Get translational energy from a set of velocities

    Parameters
    ----------
    velocities : list
        list of molecular velocities

    masses : list
        masses of the corresponding atoms

    Returns
    -------
    translational_energy_array : list
        translational energy for each of the velocities (units depending on inputs)

    """
    translational_energy_array = []
    # loop through the list of velocities
    for s in range(0, len(velocities)):
        # calculate energy # E = 0.5mv^2
        energy = 0.5 * np.sum(masses) * (velocities[s]) ** 2
        translational_energy_array.append(energy)

    return translational_energy_array


def polar_coordinates(outer_atoms, center_atoms, md_file):
    """
    Calculate the orientation of a molecule in polar coordinates
    Specific to H2O but could be adjusted

    Parameters
    ----------
    outer_atoms : list
        the atoms orbiting around a central atom (or point)

    center_atoms : list
        the atom (or point) which other atoms are orbiting around

    md_file : string
        The file to read

    Returns
    -------
    theta_array : list
        polar angle theta values for every time step degrees

    phi_array : list
        azimuthal angle phi values for every time step degrees

    rho_array : list
        radial distance for every time step Angstroms

    sin_phi_array : list
        sin of phi for every time step

    """

    # get coordinates of atoms to serve as the central or outer atoms
    outer_data = get_data(outer_atoms, md_file, 'R')
    inner_data = get_data(center_atoms, md_file, 'R')

    # containers for the angles at every time step
    theta_array = []
    phi_array = []
    rho_array = []
    sin_phi_array = []

    to_ang = 0.529177249

    # loop round every timestep
    for t in range(0, len(outer_data[0]), 1):
        # set the orientation vector to be between the two hydrogen atoms
        # calculate average H position and convert to SI
        a_vx = (float(outer_data[0][t][2]) * to_ang + float(outer_data[1][t][2]) * to_ang) / 2
        a_vy = (float(outer_data[0][t][3]) * to_ang + float(outer_data[1][t][3]) * to_ang) / 2
        a_vz = (float(outer_data[0][t][4]) * to_ang + float(outer_data[1][t][4]) * to_ang) / 2
        # find average center atom position origin
        o_x = float(inner_data[0][t][2]) * to_ang
        o_y = float(inner_data[0][t][3]) * to_ang
        o_z = float(inner_data[0][t][4]) * to_ang
        # find the difference to give direction
        x_diff = a_vx - o_x
        y_diff = a_vy - o_y
        z_diff = a_vz - o_z

        # convert to spherical coordinates with notation as used in physics
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system

        # calculate the angles
        theta = math.acos(z_diff/math.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2))
        phi = np.sign(y_diff) * math.acos(x_diff/math.sqrt(x_diff ** 2 + y_diff ** 2))

        # calculate the distance
        rho = np.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)

        # convert to degrees and append the values to the lists
        theta_array.append(math.degrees(theta))
        phi_array.append(math.degrees(phi))
        sin_phi_array.append(np.sin(phi))
        rho_array.append(rho)

    return theta_array, phi_array, rho_array, sin_phi_array


def phi_continuous(phi_list):
    """
    Manages discontinuities with polar coordinates (jumping from +180 to - 180)

    Parameters
    ----------
    phi_list : list
        phi angles in degrees

    Returns
    ----------
    phi_arr2 : list
        phi angles in degrees with discontinuities removed by adding or subtraction

    phi_arr_m
        phi angles in degrees with discontinuities removed by 'mirroring' to change the sign (probably useless)

    """

    # use a numpy array to hold the phi data
    phi_arr2 = np.array(phi_list)

    # switch the sign of the first value
    phi_arr_m = list(phi_list)
    phi_arr_m[0] = phi_arr_m[0] * -1

    # for every value in the list except the first
    for i in range(1, len(phi_arr2)):

        # find the difference in phi in the before and after steps
        difference = phi_arr2[i - 1] - phi_arr2[i]
        # add or subtract degrees when crosses a boundary i.e a sudden change beyond 300 degrees
        if difference > 350:
            phi_arr2[i:] = phi_arr2[i:] + 360
        elif difference < -350:
            phi_arr2[i:] = phi_arr2[i:] - 360

        # mirror the value by switching sign
        if phi_arr_m[i] < 0:
            phi_arr_m[i] = phi_arr_m[i] * -1

    return phi_arr2, phi_arr_m


def get_relative_orientation(outer_atoms, center_atoms, md_file):
    """
    Calculates the relative orientation of molecules by determining the bond vectors
    and dipole moments of a central atom (e.g., oxygen) with two outer atoms (e.g., hydrogen)
    over multiple time steps in a molecular dynamics (MD) simulation.

    Parameters
    ----------
    outer_atoms : list
        List of outer atoms' indices (e.g., the hydrogen atoms in water).
    center_atoms : list
        List of central atoms' indices (e.g., the oxygen atom in water).
    md_file : str
        Path to the molecular dynamics (MD) file containing positional data.

    Returns
    -------
    v1as : list of lists
        A list of bond vectors from the central atom to the first outer atom (hydrogen 1)
        for each time step. Each sublist contains the x, y, z components of the vector.
    v2as : list of lists
        A list of bond vectors from the central atom to the second outer atom (hydrogen 2)
        for each time step. Each sublist contains the x, y, z components of the vector.
    dipoles : list of lists
        A list of dipole moments (average of the two bond vectors) for each time step.
        Each sublist contains the x, y, z components of the dipole vector.
    """

    # Fetch the position data (coordinates) for the outer (e.g., hydrogen) and center (e.g., oxygen) atoms
    # from the molecular dynamics (MD) file for each time step.
    outer_data = get_data(outer_atoms, md_file, 'R')  # 'R' specifies positional data
    inner_data = get_data(center_atoms, md_file, 'R')

    # Initialize empty lists to store vectors and dipoles
    v1as = []  # List to store the first vector representing a bond (hydrogen 1 to oxygen)
    v2as = []  # List to store the second vector representing a bond (hydrogen 2 to oxygen)
    dipoles = []  # List to store the dipole moment vector

    # Loop over each time step in the data
    for t in range(0, len(outer_data[0]), 1):
        # Extract and define the 3D coordinates of the atoms at time step t:
        # 'inner_data[0]' corresponds to the oxygen atom, 'outer_data[0]' and 'outer_data[1]'
        # correspond to the two hydrogen atoms.

        # Get the coordinates of the center atom (oxygen) at time step t
        oxygen = np.array([float(inner_data[0][t][2]), float(inner_data[0][t][3]), float(inner_data[0][t][4])])

        # Get the coordinates of the first outer atom (hydrogen1) at time step t
        hydrogen1 = np.array([float(outer_data[0][t][2]), float(outer_data[0][t][3]), float(outer_data[0][t][4])])

        # Get the coordinates of the second outer atom (hydrogen2) at time step t
        hydrogen2 = np.array([float(outer_data[1][t][2]), float(outer_data[1][t][3]), float(outer_data[1][t][4])])

        # Calculate the vectors representing the bonds between oxygen and the two hydrogen's
        v1a = hydrogen1 - oxygen  # Vector from oxygen to first hydrogen
        v2a = hydrogen2 - oxygen  # Vector from oxygen to second hydrogen

        # Calculate the dipole moment by averaging the two bond vectors.
        # This represents the orientation of the water molecule based on the two hydrogen bonds.
        dipole = (v1a + v2a) / 2

        # Append the vectors and dipole to their respective lists
        v1as.append(v1a.tolist())  # Store bond vector 1 (hydrogen1 to oxygen)
        v2as.append(v2a.tolist())  # Store bond vector 2 (hydrogen2 to oxygen)
        dipoles.append(dipole.tolist())  # Store dipole moment

    # Return the lists of bond vectors and dipoles
    return v1as, v2as, dipoles


def trajectories(atoms, md_file):
    """
    Find the trajectories of a list of atoms

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    Returns
    -------
    co_xs : list of arrays
        contains the x coordinates for every time step for each atom specified

    co_ys : list of arrays
        contains the y coordinates for every time step for each atom specified

    co_zs : list of arrays
        contains the z coordinates for every time step for each atom specified

    atom_id : list
        atom IDs

    """
    # get position data from file
    data = get_data(atoms, md_file, 'R')

    to_ang = 0.529177249

    # find xyz positions of atoms of Atom list and convert bohr to angstroms
    co_xs = []
    co_ys = []
    co_zs = []
    atom_id = []
    # loop through every atom in the data
    for atom in data:
        co_x = []
        co_y = []
        co_z = []
        # append the ID of each atom
        atom_id.append(atom[0][7])
        # loop through every timestep
        for time in range(0, len(data[0]), 1):
            # convert to angstroms then append the data to lists
            co_x.append(float(atom[time][2]) * to_ang)
            co_y.append(float(atom[time][3]) * to_ang)
            co_z.append(float(atom[time][4]) * to_ang)
        # append the coordinates for each atom at every timestep
        co_xs.append(co_x)
        co_ys.append(co_y)
        co_zs.append(co_z)

    # co_x[atoms][x coordinate] etc ... atom_id[id]
    return co_xs, co_ys, co_zs, atom_id


def get_angular_velocity_direction(outer_atoms, center_atoms, md_file, atoms, masses):
    """
    Calculate the angular velocity with respect to the direction of motion
    Needs to be amended for dynamic 2D and 3D use

    Parameters
    ----------
    outer_atoms : list
        the atoms orbiting around a central atom (or point)

    center_atoms : list
        the atom (or point) which other atoms are orbiting around

    md_file : string
        The file to read

    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------

    xy_angular_velocity_max : list
        maximum angular velocity in the xy plane

    zy_angular_velocity_max : list
        maximum Angular Velocity perpendicular to direction

    zx_angular_velocity_max : list
        maximum Angular Velocity around direction of motion

    angular_velocities_atom : array of lists
        angular velocity for individual atoms

    lead : array of lists
        what atom is responsible for the maximum angular velocity result
    """

    # get position data for outer and inner atoms
    outer_data = get_data(outer_atoms, md_file, 'R')
    center_data = get_data(center_atoms, md_file, 'R')

    # define the x-axis vector
    x_axis = [1, 0, 0]

    # find the direction of motion for each step based on the center mass of the molecule
    x_dir, y_dir, z_dir = get_step_direction(atoms, md_file, masses)

    # containers for velocities in different plains
    xy_angular_velocity_max = []
    zy_angular_velocity_max = []
    zx_angular_velocity_max = []
    angular_velocities_atom = []

    # containers for leading atoms in differnt plains
    lead_counter_xy = []
    lead_counter_zy = []
    lead_counter_zx = []

    # loop through every timestep
    for t in range(0, len(x_dir), 1):
        # container for the angular velocity at a particular timestep
        angular_velocity_atom = []

        # move x-axis to be the direction of motion
        t_direction_vector = [x_dir[t], y_dir[t], 0]  # z_dir[t]] replace 0 for 3D
        rotation_matrix_xt = rotation_matrix_from_vectors(x_axis, t_direction_vector)

        # find the rotation speed for each atom
        for a in range(0, len(outer_atoms), 1):
            # recalculate points based on the new vector
            xyz_center_t = [float(center_data[0][t][2]) * 0.529177249,
                            float(center_data[0][t][3]) * 0.529177249,
                            float(center_data[0][t][4]) * 0.529177249]
            xyz_outer_t = [float(outer_data[a][t][2]) * 0.529177249,
                           float(outer_data[a][t][3]) * 0.529177249,
                           float(outer_data[a][t][4]) * 0.529177249]

            xyz_center_new_t = np.dot(rotation_matrix_xt, xyz_center_t)
            xyz_outer_new_t = np.dot(rotation_matrix_xt, xyz_outer_t)
            xyz_center_t1 = [float(center_data[0][t + 1][2]) * 0.529177249,
                             float(center_data[0][t + 1][3]) * 0.529177249,
                             float(center_data[0][t + 1][4]) * 0.529177249]
            xyz_outer_t1 = [float(outer_data[a][t + 1][2]) * 0.529177249,
                            float(outer_data[a][t + 1][3]) * 0.529177249,
                            float(outer_data[a][t + 1][4]) * 0.529177249]
            xyz_center_new_t1 = np.dot(rotation_matrix_xt, xyz_center_t1)
            xyz_outer_new_t1 = np.dot(rotation_matrix_xt, xyz_outer_t1)

            # construct and find the length of the 3 sided triangle for each axis
            center_vdis_t0 = []
            center_vdis_t1 = []
            outer_vdis_t0t1 = []
            for v in range(0, 3, 1):
                average_center = (xyz_center_new_t[v] + xyz_center_new_t1[v]) / 2
                center_vdis_t0.append(average_center - xyz_outer_new_t[v])
                center_vdis_t1.append(average_center - xyz_outer_new_t1[v])
                outer_vdis_t0t1.append(xyz_outer_new_t[v] - xyz_outer_new_t1[v])

            center_dis_t0 = []
            center_dis_t1 = []
            outer_dis_t0t1 = []
            # find distances for each plane
            for v in range(0, len(center_vdis_t0), 1):
                center_dis_t0.append(math.sqrt(center_vdis_t0[v] ** 2 + center_vdis_t0[(v + 1) % 3] ** 2))
                center_dis_t1.append(math.sqrt(center_vdis_t1[v] ** 2 + center_vdis_t1[(v + 1) % 3] ** 2))
                outer_dis_t0t1.append(math.sqrt(outer_vdis_t0t1[v] ** 2 + outer_vdis_t0t1[(v + 1) % 3] ** 2))

            angular_velocities = []
            # find arc lengths for each axis from average length and differences between vectors.
            for v in range(0, len(center_vdis_t0), 1):
                average_r = (center_dis_t0[v] + center_dis_t1[v]) / 2

                v1 = [center_vdis_t0[v], center_vdis_t0[(v + 1) % 3]]
                v2 = [center_vdis_t1[v], center_vdis_t1[(v + 1) % 3]]

                angle = np.arctan2(np.cross(v1, v2), np.dot(v1, v2))
                arc_length = angle * average_r
                angular_velocity = arc_length / 1  # distance /time
                # group all angular velocities for each time step together [xy,yz,zx]
                angular_velocities.append(angular_velocity)

            angular_velocity_atom.append(angular_velocities)
        angular_velocities_atom.append(angular_velocity_atom)

        # compare the angular velocities of each atom
        xy_max = max(angular_velocity_atom[0][0], angular_velocity_atom[1][0])
        yz_max = max(angular_velocity_atom[0][1], angular_velocity_atom[1][1])
        zx_max = max(angular_velocity_atom[0][2], angular_velocity_atom[1][2])
        xy_angular_velocity_max.append(xy_max)
        zy_angular_velocity_max.append(yz_max)
        zx_angular_velocity_max.append(zx_max)

        # compare to see if the atom leading the rotation has changed
        if angular_velocity_atom[0][0] > angular_velocity_atom[1][0]:
            lead_counter_xy.append(outer_atoms[0])
        else:
            lead_counter_xy.append(outer_atoms[1])

        if angular_velocity_atom[0][1] > angular_velocity_atom[1][1]:
            lead_counter_zy.append(outer_atoms[0])
        else:
            lead_counter_zy.append(angular_velocity_atom[1][1])

        if angular_velocity_atom[0][2] > angular_velocity_atom[1][2]:
            lead_counter_zx.append(outer_atoms[0])
        else:
            lead_counter_zx.append(outer_atoms[1])

    lead = [lead_counter_xy, lead_counter_zy, lead_counter_zx]

    return xy_angular_velocity_max, zy_angular_velocity_max, zx_angular_velocity_max, angular_velocities_atom, lead


def change_direction(atoms, md_file, masses, t_range):
    """
    determine on what time steps an atom or molecule change direction
    with respect to the x y z coordinates

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    masses : list
        relative masses of the corresponding atoms

    t_range : float
        arbitrary number of time steps around the direction change time step

    Returns
    -------
    changes_array_xyz : list
        time steps around when a change in direction event occurs relative to direction of motion

    changes_count_xyz : integer
        total number of change direction events relative to direction of motion

    around_change_co_x : list
        time steps around when a change in direction event occurs along the x-axis

    around_change_co_y : list
        time steps around when a change in direction event occurs along the y-axis

    around_change_co_z : list
        time steps around when a change in direction event occurs along the z-axis

    times_list : list of arrays
        times when the there is a direction change in x y or z (without time steps around)

    """

    # get the direction of motion for the molecule at each time step and put in a list
    x_vector_array, y_vector_array, z_vector_array = get_step_direction(atoms, md_file, masses)
    vector_arrays = [x_vector_array, y_vector_array, z_vector_array]

    # containers for:
    changes_array_xyz = []  # timestep of changes
    changes_count_xyz = []  # number of changes
    times_list = []  # number of timesteps
    # for every vector in each x y z vector list
    for v in range(0, len(vector_arrays)):
        # container to hold timestep and number of changes 1 for moving in +ve direction -1 -ve direction)
        changes_array = []
        for i in range(0, len(x_vector_array) - 1, 1):
            # compare to last step and see if it is moving in +ve or -ve directions
            if vector_arrays[v][i] > 0:
                changes_array.append(1)
            elif vector_arrays[v][i] < 0:
                changes_array.append(-1)

    # counter for number of changes
    changes_count = 0
    # for every timestep (excluding the last)
    for j in range(1, len(changes_array), 1):
        # calculate number of changes
        # if the value is not the same as the previous value then it has changed direction
        if changes_array[j] != changes_array[j - 1]:
            # increase the count by 1
            changes_count = changes_count + 1
            # record the index of when changes take place
            times_list.append(j)

    # apend the number of changes and the indexes of the change
    changes_count_xyz.append(changes_count)
    changes_array_xyz.append(changes_array)

    # to get the positions for the time frame around direction change
    # get the positions of atoms in the file
    co_xs, co_ys, co_zs, atom_id = trajectories(atoms, md_file)

    # containers for coordinates around when the direction change took place for every atom
    around_change_co_x = []
    around_change_co_y = []
    around_change_co_z = []
    # for every atom selected
    for m in range(0, len(co_xs)):
        # containers for coordinates around when the direction change took place
        around_change_x = []
        around_change_y = []
        around_change_z = []
        # for every index in the list of when the direction changes took place
        for x in range(0, len(times_list)):
            # define a time range around change event
            l_bound = times_list[x] - t_range
            u_bound = times_list[x] + t_range
            # get xyz coordinate for each atom based on these time ranges
            for n in range(int(l_bound), int(u_bound)):
                # ensure bounds are valid
                # no negative timesteps or beyond the length of the simulation
                if l_bound >= 0 and u_bound <= len(co_xs[0]):
                    # append the coordinates around when the changes toot please
                    around_change_x.append(co_xs[m][n])
                    around_change_y.append(co_ys[m][n])
                    around_change_z.append(co_zs[m][n])
        # append to the total list of coordinates
        around_change_co_x.append(around_change_x)
        around_change_co_y.append(around_change_y)
        around_change_co_z.append(around_change_z)

    return changes_array_xyz, changes_count_xyz, around_change_co_x, around_change_co_y, around_change_co_z, times_list


def rotation_matrix_from_vectors(vec1, vec2):
    """
    Finds the rotation matrix that aligns vec1 to vec2

    Parameters
    ----------
    vec1 : list
        starting vector
    vec2 : list
        new vector

    Returns
    -------
    matrix : list of lists
        when applies to vec1 aligns with vec2
    """

    # normalise the vectors a and b to ensure x y z
    a = (vec1 / np.linalg.norm(vec1)).reshape(3)
    b = (vec2 / np.linalg.norm(vec2)).reshape(3)

    # calculate the cross product to give a perpendicular angle v
    v = np.cross(a, b)

    # calculate the dot product to give alignment of vectors c
    c = np.dot(a, b)

    # calculate the sin of the angle s from the cross product of a and b
    s = np.linalg.norm(v)

    # generate the Skew-Symmetric Matrix kmat
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    # apply Rodrigues' rotation formula to give the rotation matrix
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix


def expand_points(x_points, y_points, name_points, pos_x, neg_x, pos_y, neg_y, lattice_a, lattice_b, angle):
    """
    Expands a set of points based on lattice parameters and the desired number of repetitions
    in both x and y directions, creating a periodic arrangement.

    Parameters
    ----------
    x_points : list
        List of initial x coordinates of the points.
    y_points : list
        List of initial y coordinates of the points.
    name_points : list
        List of names (or labels) associated with each point.
    pos_x : int
        Positive boundary for replication along the x-axis.
    neg_x : int
        Negative boundary for replication along the x-axis.
    pos_y : int
        Positive boundary for replication along the y-axis.
    neg_y : int
        Negative boundary for replication along the y-axis.
    lattice_a : float
        Length of the lattice along the x-axis.
    lattice_b : float
        Length of the lattice along the y-axis.
    angle : float
        Angle between the lattice vectors (in degrees).

    Returns
    -------
    new_x_points : list
        List of new x coordinates for the expanded points.
    new_y_points : list
        List of new y coordinates for the expanded points.
    new_point_names : list
        List of names for the expanded points.
    scale : float
        Scaling factor for plotting or visualization purposes.
    """

    # Calculate the acute angle in radians (adjusting from 90 degrees)
    acute_angle = np.radians(angle - 90)

    # Calculate the y-axis shift and additional x-shift based on the lattice vectors and angle
    y_axis = np.cos(acute_angle) * lattice_b  # Projection of lattice_b on y-axis
    additional_x = -np.sin(acute_angle) * lattice_b  # Component of lattice_b along the x-axis

    # Repeat the original x, y coordinates and names for the number of repetitions in the x and y directions
    x_coords = np.tile(x_points, (pos_x - neg_x, pos_y - neg_y))
    y_coords = np.tile(y_points, (pos_x - neg_x, pos_y - neg_y))
    name_coords = np.tile(name_points, (pos_x - neg_x, pos_y - neg_y))

    # Generate the x and y shifts for each repetition based on the lattice dimensions
    x_shift = np.repeat(np.arange(neg_x, pos_x), len(x_points))
    y_shift = np.repeat(np.arange(neg_y, pos_y), len(x_points))

    # Calculate the new x and y coordinates by adding lattice displacements
    new_x_points = x_coords + x_shift * lattice_a + y_shift * additional_x
    new_y_points = y_coords + y_shift * y_axis
    new_point_names = name_coords

    # Convert the new coordinates and point names from NumPy arrays to lists
    new_x_points = new_x_points.tolist()
    new_y_points = new_y_points.tolist()
    new_point_names = new_point_names.tolist()

    # Calculate a scaling factor for visualizing the grid, ensuring the boundary scale is not zero
    boundary_scale = max(pos_x - neg_x, pos_y - neg_y)
    scale = 20 / boundary_scale if boundary_scale != 0 else 1  # Avoid division by zero

    # Return the first elements of the lists along with the scale factor useful for plotting a dynamic number of atoms
    return new_x_points[0], new_y_points[0], new_point_names[0], scale


def site_analysis(atoms, md_file, lattice_a, lattice_b, angle, surface_atoms, new_point_dis, new_point_name,
                  time_array):
    """
    Evaluates what surface site an atom is above
    Note some values are hard coded in and may need to be revised depending on system
    I plan to make it more dynamic at some point based using the unit cell parameters function
    I found this function to take some time so would recommend piping the output to a text file

    Parameters
    ----------
    atoms : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    lattice_a : float
        lattice vector b length

    lattice_b : float
        lattice vector b length

    angle : float
        angle between lattice vectors a and b

    surface_atoms : list
        Atoms that make up the surface
        Can select atoms using the atom ID as appears in the .md file e.g. ["C1", "C2", "C3"]

    new_point_dis : list of lists
        distance of new points from existing point
        each item in the list consists of an x and y distance

    new_point_name : list
        corresponding identifiers for the new points

    time_array : list
        simulation time steps to include

    Returns
    -------
    s_atom_distances_array : list
        xy distance between the surface atom site and the atom

    closet_site_array : list
        list of closest site to the atom at each time step

    """

    tx, ty, tz, ta = trajectories(atoms, md_file)
    s_atom_data = get_data(surface_atoms, md_file, 'R')

    s_atom_distances_array = []
    closet_site_array = []

    acute_angle = float(angle) - 90
    angle_a = (math.cos(math.radians(acute_angle)))

    # for each time step in trajectory find out where the surface atoms are and expand
    for e in range(0, len(time_array), 1):
        t = time_array[e]

        # get atom positions and additional points at time t and put in the same arrays
        # take care with the number of surface atoms to prevent duplication
        adx, ady, names = additional_points(surface_atoms, t, md_file, new_point_dis, new_point_name)

        for k in range(0, len(surface_atoms), 1):
            adx.append(float(s_atom_data[k][t][2]) * 0.529177249)
            ady.append(float(s_atom_data[k][t][3]) * 0.529177249)
            names.append(s_atom_data[k][t][7])
        # find distance from averages
        x_dis_mean_to_point = float(tx[0][t]) - statistics.mean(adx)
        y_dis_mean_to_point = float(ty[0][t]) - statistics.mean(ady)

        # account for angle between axis
        y_axis_height = angle_a * lattice_b

        # find ratioes of distacnes
        x_dis_ratio = x_dis_mean_to_point / lattice_a
        y_dis_ratio = y_dis_mean_to_point / y_axis_height
        x_dis_ratio_rounded = round(x_dis_ratio)
        y_dis_ratio_rounded = round(y_dis_ratio)

        # find appropriate correction and round as expand points requires an interger value
        x_shift_ratio = y_dis_ratio # ratio of x and y shifts may need to change depending on system
        x_shift_correction = round(x_shift_ratio)

        pos_y = 2 + 1 * y_dis_ratio_rounded
        pos_x = 2 + x_dis_ratio_rounded + x_shift_correction
        neg_x = pos_x - 3
        neg_y = pos_y - 3

        x_ps, y_ps, p_ns, s = expand_points(adx, ady, names, pos_x, neg_x, pos_y, neg_y, lattice_a, lattice_b, angle)

        # loop through expanded points to find the site that is closest for each time
        s_atom_distances = []
        for i in range(0, len(x_ps)):
            xs_distance = tx[0][t] - x_ps[i]
            ys_distance = ty[0][t] - y_ps[i]
            s_atom_distances.append(math.sqrt(xs_distance ** 2 + ys_distance ** 2))

        min_s_atom_distance = min(s_atom_distances)
        s_atom_distances_array.append(min_s_atom_distance)

        # by creating another list of the same length as the positions can then use the index to find the site name
        min_s_atom_index = s_atom_distances.index(min_s_atom_distance)

        min_s_atom_name = p_ns[min_s_atom_index]
        closet_site_array.append(min_s_atom_name)

    return s_atom_distances_array, closet_site_array


def calculate_vacf(velocities, interval):
    """
    Calculate the velocity autocorrelation function (VACF) for a given set of velocities.

    Parameters
    ----------
    velocities : array-like
        A list or array of velocity vectors at different time steps.
    interval : int
        Time interval between samples (used to calculate the VACF at different time lags).

    Returns
    -------
    vacf : numpy array
        The normalized velocity autocorrelation function (VACF) as a function of time lag.
    """

    # Total number of time steps in the velocity data
    total_steps = len(velocities)

    # Number of intervals (or time lags) for which VACF will be calculated
    num_steps = int(total_steps / interval)

    # Initialize VACF array with zeros, size determined by the number of time steps
    vacf = np.zeros(num_steps)

    # Calculate the VACF for each time lag (t)
    for t in range(num_steps):
        autocorr = 0  # Variable to store the sum of velocity dot products
        count = 0  # Counter to track the number of contributions for normalization

        # Loop over the time origins (t0), ensuring that the time lag (t * interval) fits within the total time
        for t0 in range(total_steps - t * interval):
            v0 = velocities[t0]  # Velocity at time t0
            vt = velocities[t0 + t * interval]  # Velocity at time t0 + t * interval (the lagged time)

            # Calculate the dot product between the velocities at t0 and t0 + t * interval
            autocorr += np.dot(v0, vt)  # Accumulating the dot products for the autocorrelation
            count += 1  # Increment the number of contributions

        # Average the autocorrelation by dividing by the number of terms
        vacf[t] = autocorr / count

    # Normalize the VACF by the value at t = 0 (initial autocorrelation)
    vacf /= vacf[0]

    # Return the normalized VACF array
    return vacf


def get_all_atoms(md_file, quantity):
    """
    get a quantity of all the atoms in a md file
    (without having to specify the atom IDs)

    Parameters
    ----------
    md_file : string
        .md file of interest

    quantity  : string
        Quantity of interest "<-- R" position "<-- F" force "<-- V" velocity

    Returns
    -------
    co_x : list
        List of x quantities
    co_y : list
        List of y quantities
    co_z : list
        List of z quantities
    atom_id : list
        ID for each atom

    """

    # get coordinates of all atoms withing a specific timestep range
    # filename, number, number '<-- R' etc..

    # read the file
    with open(md_file, "r") as f:
        # container for list of quantities for all timestamps
        list0 = []

        # loop over every line in file
        for line in f:
            # if the line contains the quantity of interest
            if quantity in line:
                # split the line
                line = line.split()
                # create and ID from the atomic symbol and atom number
                identifier = line[0] + line[1]
                # append the line with the ID
                line.append(identifier)
                list0.append(line)

    # contains for xyz values with id
    co_x = []
    co_y = []
    co_z = []
    atom_id = []
    # loop through all selected atoms
    for atom in list0:
        # append the id
        atom_id.append(atom[7])
        # append the x quantity
        co_x.append(float(atom[2]))
        # append the y quantity
        co_y.append(float(atom[3]))
        # append the z quantity
        co_z.append(float(atom[4]))

    # co_x[atoms][x coordinate] etc ... atom_id[id]
    return co_x, co_y, co_z, atom_id


def calculate_msd(positions, interval=1):
    """
    Calculate the Mean Squared Displacement (MSD) and RMSD from a NumPy array of positions.

    Parameters
    ----------
    positions : np.array
        A 2D NumPy array where each row contains (x, y, z) coordinates at a given time step in Angstroms
    interval : int, optional
        The interval for calculating the MSD. Default is 1

    Returns
    -------
    msd : list
        A list where each element is the MSD for a given time lag (tau) m^2
    diffusion_coeffs_msd : list
        A list of local diffusion coefficients based on the MSD. in m^2/s
    rmsd : list
        A list where each element is the RMSD (square root of MSD)
    diffusion_coeffs_rmsd : list
        A list of local diffusion coefficients based on the RMSD
    """

    # Ensure positions are a NumPy array
    # 2D array of positions at every timestep [x1_t0, y1_t0, z1_t0], [x2_t0, y2_t0, z2_t0]
    positions = np.array(positions) * 1E-10  # correct to meters

    # Lists to store MSD and RMSD values for every time step
    msd = []
    rmsd = []

    # Calculate MSD for each time lag (tau)
    # for every timestep except the first (with interval control for interest)
    for tau in range(1, len(positions), interval):
        # calculate the displacements
        # by subtracting the position at every time form the position at every time + tau
        displacements = positions[tau:] - positions[:-tau]
        # square the displacements
        squared_displacements = displacements ** 2
        # sum over x y and z coordinates for each time step
        summed_squared_displacements = np.sum(squared_displacements, axis=1)
        # calculate the mean squared displacement
        mean_squared_displacement = np.mean(summed_squared_displacements)
        # append the msd
        msd.append(mean_squared_displacement)
        # append the rmsd
        rmsd.append(math.sqrt(mean_squared_displacement))

    # Calculate local slopes for the diffusion coefficient
    diffusion_coeffs_msd = []
    diffusion_coeffs_rmsd = []

    for i in range(1, len(msd)):
        # MSD-based diffusion coefficient
        delta_msd = msd[i] - msd[i - interval]
        delta_time = interval * 1E-15  # Time interval in seconds
        local_slope = delta_msd / delta_time
        local_diffusion_coeff = local_slope / 6  # 6 for 3D diffusion
        diffusion_coeffs_msd.append(local_diffusion_coeff)

        # RMSD-based diffusion coefficient
        delta_rmsd = rmsd[i] - rmsd[i - interval]
        local_slope_rmsd = delta_rmsd / delta_time
        local_diffusion_coeff_rmsd = local_slope_rmsd / 6  # 6 for 3D diffusion
        diffusion_coeffs_rmsd.append(local_diffusion_coeff_rmsd)

    return msd, diffusion_coeffs_msd, rmsd, diffusion_coeffs_rmsd


def calculate_tamsd(positions, time_lags):
    """
    Calculate the Time-Averaged Mean Squared Displacement (TAMSD) for 3D positions from a list.
    Reduces noise in the msd by considering total time of the observation

    Parameters:
    positions (list of lists): List of molecule positions at different timesteps (N_timesteps, 3).
    time_lags (ndarray): Array of time lags (Δt) at which to calculate TAMSD.

    Returns:
    tamsd (ndarray): Array of TAMSD values for each time lag.
    """

    # Ensure positions are a NumPy array
    # 2D array of positions at every timestep [x1_t0, y1_t0, z1_t0], [x2_t0, y2_t0, z2_t0]
    positions = np.array(positions) * 1E-10  # Convert list of lists to NumPy array
    # get the number of timestamps and dimensions
    n_timesteps, n_dimensions = positions.shape
    # initialise the tamsd array with 0 for each of the time lags
    tamsd = np.zeros(len(time_lags))

    # Convert time lags from femto-seconds to seconds
    time_lags_in_seconds = np.array(time_lags) * 1E-15  # Convert fs to seconds

    # Loop over all time lags
    for i, delta_t in enumerate(time_lags):
        squared_displacements = []

        # for every time in the time lags
        for t in range(n_timesteps - delta_t):
            # Calculate squared displacement for each time lag
            displacement = positions[t + delta_t] - positions[t]
            squared_displacements.append(np.sum(displacement ** 2))  # Sum of x^2 + y^2 + z^2

        # Time-average the squared displacements
        tamsd[i] = np.mean(squared_displacements)

    # Remove the first value as it should be 0
    diffusion_coefficients = tamsd[1:] / (6 * time_lags_in_seconds[1:])

    return tamsd, diffusion_coefficients


def effective_diffusion_constant_msd(msd, times, alpha):
    """
    Calculate the effective diffusion constant based on msd times and alpha.
    Essentially fit a curve to the msd to account for supper and sub diffusion

    Parameters:
    msd : list
        mean squared displacement Angstroms^2

    times : list
        timesteps of the simulation

    alpha : float
        The scaling exponent from the MSD power-law fit.

    Returns:
    np.ndarray: Effective diffusion constants at the given timesteps.
    """

    # Adjust based on the exponent alpha
    times = np.asarray(times[1:])
    # Need to remove the first step as =0 by definition
    msd = np.asarray(msd[1:])

    if alpha == 1:
        # Simple diffusion (standard case)
        return msd / (6 * times)
    elif alpha < 1:
        # Sub-diffusion case
        return msd / (6 * (times) ** alpha)
    else:
        # Super-diffusion case (generalized approach, can be more complex)
        return msd / (6 * (times) ** alpha)


def effective_diffusion_constant_C_fit(C_fit, times, alpha):
    """
    Calculate the effective diffusion constant based on C_fit, time, and alpha
    (opposed to the msd)

    Parameters:
    C_fit (float): The pre-factor from the MSD power-law fit
    times (np.ndarray): Array of time points.
    alpha (float): The scaling exponent from the MSD power-law fit

    Returns:
    np.ndarray: Effective diffusion constants at the given time points
    """
    # Convert times to a numpy array if it's not already
    times = np.asarray(times)

    # Adjust based on the exponent alpha
    if alpha == 1:
        # Normal diffusion (standard case)
        D_eff = C_fit / 6  # Constant diffusion
    else:
        # Anomalous diffusion (sub-diffusion or super-diffusion)
        D_eff = (C_fit * alpha * (times ** (alpha - 1))) / 6

    return D_eff


def calculate_heat_capacity_system(md_file):
    """
    Calculate the heat capacity of the system over time.

    Parameters
    ----------
    md_file : string
        The file to read the data from.

    Returns
    -------
    times : array-like
        Array of time steps.
    heat_capacities : array-like
        Heat capacities at each time step.
    """
    # Get energy and temperature data
    energy = get_energy(md_file)
    temperature = get_temperature(md_file)

    # includes all 3 forms of energy total, Hamiltonian and kinetic energy
    energy = energy[3]
    temperature = temperature

    num_steps = len(energy)
    heat_capacities = []

    for i in range(num_steps):
        # Calculate the mean energy up to the current time step
        avg_energy = np.mean(energy[:i + 1])

        # Calculate energy fluctuations up to the current time step
        energy_fluctuation = np.mean((energy[:i + 1] - avg_energy) ** 2)

        # Calculate the temperature at the current time step
        temp = temperature[i]

        # Calculate the heat capacity at the current time step
        if temp != 0:
            heat_capacity = energy_fluctuation / (temp ** 2)
            heat_capacities.append(heat_capacity)
        else:
            heat_capacities.append(np.nan)  # Handle zero temperature case

    # Create an array of time steps
    times = np.arange(num_steps)  # Adjust if you have specific time intervals

    return times, np.array(heat_capacities)


def calculate_heat_capacity_molecule(atoms, md_file, masses):
    """
    Calculate the heat capacity of the system over time for a molecule based on the atomic velocities and temperature.

    Parameters
    ----------
    atoms : list
        List of atoms to consider when calculating the heat capacity.
        The atoms can be selected using the atom IDs as returned by the 'get_data' function, e.g., ["H1", "H2", "H3"].

    md_file : string
        The file from which molecular dynamics (MD) data is read (likely containing positions, velocities, etc.).

    masses : list
        A list of masses corresponding to the atoms in the same order as provided in the `atoms` list.

    Returns
    -------
    times : array-like
        An array of time steps for the simulation.

    heat_capacities : array-like
        The calculated heat capacity at each time step.
    """

    # Get the velocity data (x_co, y_co, z_co) for the atoms of interest
    # `center_of_velocity` likely calculates the velocity of the system's center of mass over time
    x_co, y_co, z_co = center_of_velocity(atoms, md_file, masses)

    # Convert atomic units of velocity to meters per second (SI units)
    au_to_m_per_s = 2.18769126364e6  # Conversion factor from atomic units to m/s

    # List to store the kinetic energy at each time step
    kinetic_energies = []

    # Sum the masses to get the total mass of the system
    mass = sum(masses)

    # Loop through the velocities in the x, y, and z directions
    for vx, vy, vz in zip(x_co, y_co, z_co):
        velocity = [vx, vy, vz]

        # Convert velocity to a NumPy array and apply the conversion factor to m/s
        velocity = np.array(velocity) * au_to_m_per_s

        # Calculate the magnitude (norm) of the velocity vector
        v_magnitude = np.linalg.norm(velocity)

        # Calculate the kinetic energy for this time step: (1/2) * m * v^2
        kinetic_energy = 0.5 * mass * v_magnitude ** 2

        # Append the calculated kinetic energy to the list
        kinetic_energies.append(kinetic_energy)

    # `energy` stores the total kinetic energy of the system at each time step
    energy = kinetic_energies

    # Get temperature data from the MD file (assumed to be pre-calculated or stored in the file)
    temperature = get_temperature(md_file)

    # The number of time steps in the temperature data
    num_steps = len(temperature)

    # Initialize a list to store heat capacity values at each time step
    heat_capacities = []

    # Loop over each time step to calculate the heat capacity
    for i in range(num_steps):
        # Calculate the average energy up to the current time step
        avg_energy = np.mean(energy[:i + 1])

        # Calculate the energy fluctuation (variance in energy)
        energy_fluctuation = np.mean((energy[:i + 1] - avg_energy) ** 2)

        # Get the temperature at the current time step
        temp = temperature[i]

        # Calculate the heat capacity if temperature is non-zero
        if temp != 0:
            heat_capacity = energy_fluctuation / (temp ** 2)  # Formula: C = (delta(E)^2) / T^2
            heat_capacities.append(heat_capacity)
        else:
            # If temperature is zero, append NaN to handle this case
            heat_capacities.append(np.nan)

    # Create an array of time steps (assuming equal time intervals; adjust if needed)
    times = np.arange(num_steps)  # Or provide specific time values if available

    # Return the times and the heat capacities as NumPy arrays
    return times, np.array(heat_capacities)


def calculate_deviation(trajectory):
    """
    Calculate the deviation of a molecule's position from its trajectory in the previous timestep.

    Parameters
    ----------
    trajectory : np.ndarray
        A NumPy array of shape (num_steps, 3) representing the molecule's trajectory.
        Each row corresponds to the (x, y, z) coordinates of the molecule at a given timestep.

    Returns
    -------
    deviations : np.ndarray
        A NumPy array of shape (num_steps - 1,) representing the deviation of the molecule's
        position from the previous timestep for each timestep.
    """

    # Calculate the displacement vector between consecutive timesteps
    displacements = trajectory[1:] - trajectory[:-1]

    # Calculate the magnitude of the displacement vectors to get the deviations
    deviations = np.linalg.norm(displacements, axis=1)

    return deviations


def power_law(t, c, alpha):
    """
    for curve fitting
    """
    return c * t**alpha


def gaussian(x, a, mu, sigma):
    """
    for gaussian fitting
    """
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


def compute_free_energies(probs, temperature=150):
    """
    Compute free energies from normalized probabilities using the Boltzmann equation.

    Args:
    - probabilities (ndarray): Array of probabilities.
    - temperature (float): Temperature in Kelvin.

    Returns:
    - free_energies (ndarray): Array of computed free energies.
    """

    # Define Boltzmann constant
    kB = 8.617333262145e-5  # Boltzmann constant in eV/K

    # Normalize probabilities by their maximum value
    normalized_probabilities = np.asarray(probs) / np.max(probs)

    # Compute free energies
    free_energies = -kB * temperature * np.log(normalized_probabilities)

    return free_energies


def calculate_potential_of_mean_force(mol_atoms, md_file, num_timesteps=3000):
    """
    Compute the potential of mean force (PMF) for a molecule.

    Parameters:
    mol_atoms : list
        List of atoms in the molecule.
    md_file : str
        File containing molecular dynamics data.
    num_timesteps : int, optional
        Number of timesteps to consider. If None, all timesteps are used.

    Returns:
    - pmf: Array of computed free energies (PMF).
    - mean_forces: Array of mean force magnitudes.
    - gradient_pmf: Gradient of the PMF.
    """

    # Get force data from file, structure: data[atom][time][columns]
    data = get_data(mol_atoms, md_file, 'F')

    # Organize data in array
    data_array = np.asarray(data)

    # Convert only the force columns (2:5) from strings to floats
    forces_array = data_array[:, :, 2:5].astype(float)

    # Number of time steps
    n_time_steps = forces_array.shape[1]

    # If num_timesteps is specified, limit the number of time steps
    if num_timesteps is not None:
        n_time_steps = min(num_timesteps, n_time_steps)

    # Set up array of timesteps and convert from fs to s
    time_values = np.arange(n_time_steps) * 1E-15

    # Initialize array to store mean force magnitudes
    mean_forces = np.zeros(n_time_steps - 1)

    # Calculate mean force magnitudes for each time step
    for t in range(n_time_steps - 1):
        # Convert from Au to N (using the appropriate conversion factor)
        forces_atoms = forces_array[:, t, :] * 8.2387235038E-8

        # Calculate force magnitudes for each atom at current time step
        forces_total = np.linalg.norm(forces_atoms, axis=1)

        # Average forces over all atoms
        mean_forces[t] = np.mean(forces_total)

    # Calculate the time intervals (differences between consecutive time steps)
    delta_time = np.diff(time_values)

    # Perform numerical integration using trapezoidal rule to get PMF
    pmf = -np.cumsum(mean_forces * delta_time)

    # Compute gradient of the PMF with respect to time
    gradient_pmf = np.gradient(pmf, time_values[1:])

    return pmf, mean_forces, gradient_pmf


def calculate_distances_sum(atoms, md_file, masses):
    """
    Calculate the distance from the start for each timestep and the total distance traveled.

    Parameters:
    atoms : int
        Number of atoms in the molecule.
    md_file : str
        Path to the file containing molecular dynamics data.

    Returns:
    - distances_from_start: Array of distances from the start for each timestep.
    - cumulative_distances: Array of cumulative distances traveled over all timesteps.
    """

    # Get center of mass
    x_cos, y_cos, z_cos = center_of_mass(atoms, md_file, masses)

    # Convert coordinates to numpy arrays
    x_coords = np.asarray(x_cos)
    y_coords = np.asarray(y_cos)
    z_coords = np.asarray(z_cos)

    # Start position at time step 0
    start_position = np.array([x_coords[0], y_coords[0], z_coords[0]])

    # Calculate distances from the start for each timestep
    distances_from_start = np.sqrt((x_coords - start_position[0]) ** 2 +
                                    (y_coords - start_position[1]) ** 2 +
                                    (z_coords - start_position[2]) ** 2)

    # Calculate distances traveled between timesteps
    delta_x = np.diff(x_coords)
    delta_y = np.diff(y_coords)
    delta_z = np.diff(z_coords)
    distances_traveled = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)

    # Calculate cumulative distances traveled
    cumulative_distances = np.zeros(len(distances_traveled) + 1)  # +1 for the starting position
    cumulative_distances[1:] = np.cumsum(distances_traveled)  # Fill from index 1 onwards

    return distances_from_start, cumulative_distances


def calculate_mean_confidence_interval(data, confidence=0.98):
    g = len(data)
    mean = np.mean(data)
    stderr = stats.sem(data)  # Standard error of the mean
    interval = stats.t.interval(confidence, g - 1, loc=mean, scale=stderr)
    return mean, interval[1] - mean  # Return mean and upper bound of CI


def linear_regression(x, y):
    """
    Perform linear regression

    Parameters:
    x : list
        Independent variable data
    y : list
        Dependent variable data

    Returns:
    - slope (float): Slope of the regression line
    - intercept (float): Intercept of the regression line
    - r_squared (float): goodness of fit of the regression
    - z_score (float): standard deviations from the mean
    - p_value (float): strength of association between variables
        (< 0.05 indicates that the coefficient is statistically significant)
    """
    # perform liner regression on x and y
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    # caculated r^2
    r_squared = r_value**2
    # Calculate Z-score for the slope
    z_score = slope / std_err

    return slope, intercept, r_squared, z_score, p_value


def multivariable_regression(X, y):
    """
    Perform multivariable regression

    Parameters:
    X : array
        2D array of independent variable data (shape: [n_samples, n_features])
    y : list
        Dependent variable data

    Returns:
    - slope (float): Slope of the regression line
    - intercept (float): Intercept of the regression line
    - r_squared (float): goodness of fit of the regression
    - z_score (float): standard deviations from the mean
    - p_value (float): strength of association between variables
        (< 0.05 indicates that the coefficient is statistically significant)
    """

    # Calculate coefficients and residuals using the least squares method
    # rank mesures the number of independent rows and columns in the matrix
    # s quantifies the variance is contained in the data
    coefficients, residuals, rank, s = np.linalg.lstsq(X, y)

    # Predictions by dot product between the matrix X and coefficient vector
    y_pred = X @ coefficients

    # Calculate R-squared
    ss_total = np.sum((y - np.mean(y)) ** 2)
    ss_residual = np.sum((y - y_pred) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Calculate standard errors
    n = X.shape[0]  # number of samples
    p = X.shape[1]  # number of parameters (including intercept)
    residual_std_error = np.sqrt(ss_residual / (n - p))

    # Calculate t-statistics and p-values
    t_stats = coefficients / residual_std_error
    p_values = [2 * (1 - stats.t.cdf(np.abs(t), n - p)) for t in t_stats]

    return coefficients, coefficients[0], r_squared, p_values
