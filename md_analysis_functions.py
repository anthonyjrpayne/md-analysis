"""
Module contains functions used to analyse data from CASTEP .md files
Could be adapted to accept position data from other sources
"""

import numpy as np
import math
import statistics


# contains functions used to analyse data from castep .md files
def get_data(atom_list, md_file, value):
    """
    Reads the .md file and returns a specified set of data

    Parameters
    ----------
    atom_list : list
        Atoms of interest to get the data for
        Can select atoms using the atom ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    md_file : string
        The file to read

    value : string
        The data set to return
        'R' Position, 'V' Velocity, 'F' Force

    Returns
    -------
    data : list of arrays
        Data for each atom at every time step
        In the format data[atom][time][Atom, AtomNumber, x,y,z, value, ID]
    """

    data = []
    for atom in atom_list:
        with open(md_file, "r") as f:
            list0 = []
            for line in f:
                if value in line:
                    line = line.split()
                    identifier = line[0] + line[1]
                    if atom == identifier:
                        line.append(identifier)
                        list0.append(line)

            data.append(list0)
    return data


def additional_points(atoms, time, md_file, new_point_dis, new_point_name):
    """
    Reads the .md file generates new points in the xy plain relative to specified atoms atoms

    Parameters
    ----------
    atoms : list
        Atoms of interest selected by ID as appears in the .md file e.g. ["H1", "H2", "H3"]

    time : integer
        time step number of interest

    md_file : string
        The file to read

    new_point_dis : list of lists
        distance of new points from existing point
        each item in the list consists of an x and y distance

    new_point_name : list
        corresponding identifiers for the new points

    Returns
    -------
    new_x_point_list : list
        list of new x coordinates

    new_y_point_list : list
        list of new y coordinates

    new_point_name_list : list
        list of identifiers for each new point
    """
    new_x_point_list = []
    new_y_point_list = []
    new_point_name_list = []

    data = get_data(atoms, md_file, "R")

    for i in range(0, len(data), 1):
        for j in range(0, len(new_point_dis)):
            # taken the opportunity to convert to angstroms here by multiplying by 0.529177249
            new_x = float(data[i][time][2]) * 0.529177249 + new_point_dis[j][0]
            new_y = float(data[i][time][3]) * 0.529177249 + new_point_dis[j][1]
            new_x_point_list.append(new_x)
            new_y_point_list.append(new_y)
            new_point_name_list.append(new_point_name[j])

    return new_x_point_list, new_y_point_list, new_point_name_list


def atom_distances(md_file, relations):
    """
     Rough estimate the the internal area and volume of a 3 atom molecule

     Parameters
     ----------
     relations : array of lists
         Used to create the shape of the molecule
         Each item a list of 2 values e.g. [['H1', 'H2'], ['H1', 'O1'], ['H2', 'O1']]

     md_file : string
         The file to read

     Returns
     -------
     distances_array : array of lists
         Distances between each relation for every time step
         distances_array[relation][time][distance]

     relation_array : list
         Order of corresponding relation IDs
     """
    # get atom list from relations
    atoms = []
    for relation in relations:
        atoms.append(relation[0])
        atoms.append(relation[1])
    atom_list = list(set(atoms))

    from get_data import get_data
    data = get_data(atom_list, md_file, 'R')

    distances = []
    distances_array = []
    relation_array = []
    # First create ID for the relationship e.g. H1-H2
    for relation in relations:
        relation_id = relation[0] + '-' + relation[1]
        relation_array.append(relation_id)

        # Then loop though the data to see if the ID matches both the atoms in the relationship
        for atom1 in range(0, len(data), 1):
            if relation[0] == data[atom1][0][7]:
                for atom2 in range(0, len(data), 1):
                    if relation[1] == data[atom2][0][7]:

                        # If a match then subtract distances, calculate the 3D distance and convert to Angstroms
                        distances = []
                        for time in range(0, len(data[0]), 1):
                            x_distance = float(data[atom1][time][2]) - float(data[atom2][time][2])
                            y_distance = float(data[atom1][time][3]) - float(data[atom2][time][3])
                            z_distance = float(data[atom1][time][4]) - float(data[atom2][time][4])
                            distance = math.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
                            distances.append(distance * 0.529177249)
        distances_array.append(distances)

    # distances_array[relation][time][distance] relation_array[relation]
    return distances_array, relation_array


def find_area(relations, md_file):

    """
    Rough estimate the the internal area and volume of a 3 atom molecule

    Parameters
    ----------
    relations : array of lists
        Used to create the shape of the molecule
        Each item a list of 2 values e.g. [['H1', 'H2'], ['H1', 'O1'], ['H2', 'O1']]

    md_file : string
        The file to read

    Returns
    -------
    area_list : list
        Area in the xy plain that is occupied by the molecule

    volume_list : list
        Volume that is occupied by the molecule
    """

    atoms = []
    for relation in relations:
        atoms.append(relation[0])
        atoms.append(relation[1])
    atom_list = list(set(atoms))

    from get_data import get_data
    data = get_data(atom_list, md_file, 'R')

    distances_array = []
    relation_array = []
    # First create ID for the relationship e.g. H1-H2
    for relation in relations:
        relation_id = relation[0] + '-' + relation[1]
        relation_array.append(relation_id)

        # Then loop though the data to see if the ID matches both the atoms in the relationship
        for atom1 in range(0, len(data), 1):
            if relation[0] == data[atom1][0][7]:
                for atom2 in range(0, len(data), 1):
                    if relation[1] == data[atom2][0][7]:

                        # If a match then subtract distances, calculate the 3D distance and convert to Angstroms
                        distances = []
                        for time in range(0, len(data[0]), 1):
                            x_distance = float(data[atom1][time][2]) - float(data[atom2][time][2])
                            y_distance = float(data[atom1][time][3]) - float(data[atom2][time][3])
                            distance = math.sqrt(x_distance ** 2 + y_distance ** 2)
                            distances.append(distance * 0.529177249)
                        distances_array.append(distances)

    area_list = []
    volume_list = []
    for t in range(0, len(distances_array[0])):

        d1s = distances_array[0][t]
        d2s = distances_array[1][t]
        d3s = distances_array[2][t]

        # Heron's Formula to find the area
        ss = (d1s + d2s + d3s) / 2
        areas = math.sqrt(ss * (ss - d1s) * (ss - d2s) * (ss - d3s))
        volume_list.append(areas)

    # distances_array[relation][time][distance] relation_array[relation]
    return area_list, volume_list


def center_of_mass(atoms, md_file, masses):
    """
    gives the coordinates of the center of mass for a molecule

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
        center of mass x coordinate for every time step

    y_co : list
        center of mass y coordinate for every time step

    z_co : list
        center of mass z coordinate for every time step
    """

    data = get_data(atoms, md_file, 'R')
    # calculate total mass of molecule
    total_mass = 0
    for mass in masses:
        total_mass = total_mass + mass

    # calculate center of mass at each time step
    x_co = []
    y_co = []
    z_co = []
    for t in range(0, len(data[0]), 1):
        center_x = 0
        center_y = 0
        center_z = 0
        # summation of position * mass followed by conversion to angstroms and division of total mass
        for i in range(0, len(atoms), 1):
            center_x = center_x + (float(data[i][t][2]) * masses[i])
            center_y = center_y + (float(data[i][t][3]) * masses[i])
            center_z = center_z + (float(data[i][t][4]) * masses[i])
        x_co.append(center_x * 0.529177249 / total_mass)
        y_co.append(center_y * 0.529177249 / total_mass)
        z_co.append(center_z * 0.529177249 / total_mass)

    return x_co, y_co, z_co


def find_closest(input_list, input_value):
    """
    finds closest value in a list

    Parameters
    ----------
    input_list: array of floats
        list to search through

    input_value: float
        value to find closest of

    Returns
    -------
    output_value: float
        the closest value to the input_value that is in input_list
    """
    arr = np.asarray(input_list)
    i = (np.abs(arr - input_value)).argmin()
    output_value = arr[i]
    return arr[output_value]


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
    reference_energy = mol_energy + sub_energy

    data = []
    with open(md_file, "r") as f:
        for line in f:
            if '<-- E' in line:
                line = line.split()
                data.append(line)

    to_eV = 27.211324570273

    total_energy = []
    for i in range(0, len(data) - 10, 1):
        total_energy.append(float(data[i][0]) * to_eV)

    adsorption_energy_array = []
    # calculate adsorption energy = together - separate
    for j in range(0, len(total_energy)):
        e_ads = total_energy[j] - reference_energy
        adsorption_energy_array.append(e_ads)

    return adsorption_energy_array


def angle_between(a, b):
    """
    Calculates angle between two vectors

    Parameters
    ----------
    a : list
        vector a matrix

    b : list
        vector b matrix

    Returns
    -------
    angle:
        angle between each vector

    """
    angle = np.arctan2(np.cross(a, b), np.dot(a, b))

    return angle


def get_energy(md_file):
    """
    Calculates angle between two vectors

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

    data = []
    with open(md_file, "r") as f:
        for line in f:
            if '<-- E' in line:
                line = line.split()
                data.append(line)

    time = list(range(0, len(data) - 10))

    to_eV = 27.211324570273

    total_energy = []
    hamiltonian_energy = []
    kinetic_energy = []
    for i in range(0, len(data) - 10, 1):
        total_energy.append(float(data[i][0]) * to_eV)
        hamiltonian_energy.append(float(data[i][1]) * to_eV)
        kinetic_energy.append(float(data[i][2]) * to_eV)

    # Time, Total, Hamiltonian, Kinetic
    return time, total_energy, hamiltonian_energy, kinetic_energy


def get_furthest_distance(atoms, md_file, length):
    """
    Get total distance traveled in simulation

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

    data = get_data(atoms, md_file, "R")

    x_distance = abs(float(data[0][0][2]) - float(data[0][length][2])) * 0.529177249
    y_distance = abs(float(data[0][0][3]) - float(data[0][length][3])) * 0.529177249
    xy_distance = math.sqrt(x_distance ** 2 + y_distance ** 2)
    xy_distance_time = xy_distance

    sum_of_distances = 0
    for t in range(0, length - 1):
        x_distance1 = abs(float(data[0][t][2]) - float(data[0][t + 1][2])) * 0.529177249
        y_distance1 = abs(float(data[0][t][3]) - float(data[0][t + 1][3])) * 0.529177249
        xy_distance1 = math.sqrt(x_distance1 ** 2 + y_distance1 ** 2)
        sum_of_distances = sum_of_distances + xy_distance1

    return xy_distance_time, sum_of_distances


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
    x_vector_array : list
        value of x vector for every time step

    y_vector_array : list
        value of y vector for every time step

    y_vector_array : list
        value of z vector for every time step
    """
    x_cos, y_cos, z_cos = center_of_mass(atoms, md_file, masses)

    x_vector_array = []
    y_vector_array = []
    z_vector_array = []
    for i in range(0, len(x_cos) - 10):
        x_vector = x_cos[i + 1] - x_cos[i]
        x_vector_array.append(x_vector)
        y_vector = y_cos[i + 1] - y_cos[i]
        y_vector_array.append(y_vector)
        z_vector = z_cos[i + 1] - z_cos[i]
        z_vector_array.append(z_vector)

    return x_vector_array, y_vector_array, z_vector_array


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
    data = []
    with open(md_file, "r") as f:
        for line in f:
            if '<-- T' in line:
                line = line.split()
                data.append(line)

    temperatures = []
    for i in range(0, len(data) - 10, 1):
        temperatures.append(float(data[i][0])*3.1577464E5)

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
    atom_xy_vs_array : list of arrays
        the velocity in the xy plain for every atom for every time step
        atom_xy_vs_array[atom][time]

    atom_xyz_vs_array : list of arrays
        the velocity in the 3D for every atom for every time step
        atom_xyz_vs_array[atom][time]

    mol_xy_array : list
        the velocity in the xy plain for the molecule every time step

    mol_xyx_array : list
        the velocity in the xyz plain for the molecule every time step
    """
    # data[atom][time][Atom, AtomNumber, x,y,z]
    data = get_data(atoms, md_file, 'R')

    atom_xy_vs_array = []
    atom_xyz_vs_array = []
    mol_xy_array = []
    mol_xyz_array = []
    for a in range(0, len(data)):
        atom_xy_vs = []
        atom_xyz_vs = []
        for t in range(0, len(data[0]) - 10):
            x_velocity = (float(data[a][t][2]) - float(data[a][t + 1][2])) / 1E10
            y_velocity = (float(data[a][t][3]) - float(data[a][t + 1][3])) / 1E10
            z_velocity = (float(data[a][t][4]) - float(data[a][t + 1][4])) / 1E10

            xy_velocity = (math.sqrt(x_velocity ** 2 + y_velocity ** 2) * 0.529177249) / 1E-15
            xyz_velocity = (math.sqrt(x_velocity ** 2 + y_velocity ** 2 + z_velocity ** 2)) * 0.529177249 / 1E-15
            atom_xy_vs.append(xy_velocity)
            atom_xyz_vs.append(xyz_velocity)

        atom_xy_vs_array.append(atom_xy_vs)
        atom_xyz_vs_array.append(atom_xyz_vs)

    # x[t] y[t] and z[t]
    x_cm, y_cm, z_cm = center_of_mass(atoms, md_file, masses)

    for c in range(0, len(x_cm) - 10):
        vx_cm = (x_cm[c] - x_cm[c + 1]) / 1E10
        vy_cm = (y_cm[c] - y_cm[c + 1]) / 1E10
        vz_cm = (z_cm[c] - z_cm[c + 1]) / 1E10

        vxy_cm = (math.sqrt(vx_cm ** 2 + vy_cm ** 2)) / 1E-15
        vxyz_cm = (math.sqrt(vx_cm ** 2 + vy_cm ** 2 + vz_cm ** 2)) / 1E-15
        mol_xy_array.append(vxy_cm)
        mol_xyz_array.append(vxyz_cm)

    # atom arrays [atom][time]
    # molecule arrays [time]
    return atom_xy_vs_array, atom_xyz_vs_array, mol_xy_array, mol_xyz_array


def get_translational_energy(velocities, masses):
    """
    Get translational energy from a set of velocities

    Parameters
    ----------
    velocities : list
        list of molecular velocities

    masses : list
        relative masses of the corresponding atoms

    Returns
    -------
    translational_energy_array : list
        translational energy for each of the velocities

    """
    # E = 0.5mv^2
    translational_energy_array = []
    for s in range(0, len(velocities)):
        energy = 0.5 * np.sum(masses) * (velocities[s]) ** 2
        translational_energy_array.append(energy)

    return translational_energy_array


def polar_coordinates(outer_atoms, center_atoms, md_file):
    """
    Calculate the orientation of a molecule in polar coordinates

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
        polar angle theta values for every time step

    phi_array : list
        azimuthal angle phi values for every time step

    rho_array : list
        radial distance for every time step

    """

    outer_data = get_data(outer_atoms, md_file, 'R')
    inner_data = get_data(center_atoms, md_file, 'R')

    length = len(outer_data[0]) - 1
    theta_array = []
    phi_array = []
    rho_array = []
    # set the orientation vector to be between the two hydrogen atoms
    for t in range(0, length, 1):
        # calculate average H position
        a_vx = (float(outer_data[0][t][2]) * 0.529177249 + float(outer_data[1][t][2]) * 0.529177249) / 2
        a_vy = (float(outer_data[0][t][3]) * 0.529177249 + float(outer_data[1][t][3]) * 0.529177249) / 2
        a_vz = (float(outer_data[0][t][4]) * 0.529177249 + float(outer_data[1][t][4]) * 0.529177249) / 2
        # find average O position origin
        o_x = float(inner_data[0][t][2]) * 0.529177249
        o_y = float(inner_data[0][t][3]) * 0.529177249
        o_z = float(inner_data[0][t][4]) * 0.529177249
        # find the difference to give direction
        x_diff = a_vx - o_x
        y_diff = a_vy - o_y
        z_diff = a_vz - o_z

        # convert to spherical coordinates as used in physics
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system

        theta = math.acos(z_diff/math.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2))
        phi = np.sign(x_diff) * math.acos(y_diff/math.sqrt(x_diff ** 2 + y_diff ** 2))

        rho = np.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)
        theta_array.append(theta)
        phi_array.append(phi)
        rho_array.append(rho)

    return theta_array, phi_array, rho_array


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
    data = get_data(atoms, md_file, 'R')
    length = len(data[0]) - 2

    # find xyz positions of atoms of Atom list and convert bohr to angstroms
    co_xs = []
    co_ys = []
    co_zs = []
    atom_id = []
    for atom in data:
        co_x = []
        co_y = []
        co_z = []
        atom_id.append(atom[0][7])
        for time in range(0, length, 1):
            co_x.append(float(atom[time][2]) * 0.529177249)
            co_y.append(float(atom[time][3]) * 0.529177249)
            co_z.append(float(atom[time][4]) * 0.529177249)
        co_xs.append(co_x)
        co_ys.append(co_y)
        co_zs.append(co_z)

    # co_x[atoms][x coordinate] etc ... atom_id[id]
    return co_xs, co_ys, co_zs, atom_id


def get_angular_velocity_direction(outer_atoms, center_atoms, md_file, atoms, masses):
    """
    Calculate the angular velocity with respect to the direction of motion

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

    outer_data = get_data(outer_atoms, md_file, 'R')
    center_data = get_data(center_atoms, md_file, 'R')
    length = len(outer_data[0]) - 10

    x_axis = [1, 0, 0]

    x_dir, y_dir, z_dir = get_step_direction(atoms, md_file, masses)

    xy_angular_velocity_max = []
    zy_angular_velocity_max = []
    zx_angular_velocity_max = []
    angular_velocities_atom = []
    lead_counter_xy = []
    lead_counter_zy = []
    lead_counter_zx = []
    for t in range(0, length, 1):
        angular_velocity_atom = []

        # move x axis to be the direction of motion
        t_direction_vector = [x_dir[t], y_dir[t], 0]  # z_dir[t]]
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
            xyz_center_t1 = [float(center_data[0][t+1][2]) * 0.529177249,
                             float(center_data[0][t+1][3]) * 0.529177249,
                             float(center_data[0][t+1][4]) * 0.529177249]
            xyz_outer_t1 = [float(outer_data[a][t+1][2]) * 0.529177249,
                            float(outer_data[a][t+1][3]) * 0.529177249,
                            float(outer_data[a][t+1][4]) * 0.529177249]
            xyz_center_new_t1 = np.dot(rotation_matrix_xt, xyz_center_t1)
            xyz_outer_new_t1 = np.dot(rotation_matrix_xt, xyz_outer_t1)

            # construct and find the length of the 3 sided triangle for each axis
            center_vdis_t0 = []
            center_vdis_t1 = []
            outer_vdis_t0t1 = []
            for v in range(0, 3, 1):
                average_center = (xyz_center_new_t[v] + xyz_center_new_t1[v])/2
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
                average_r = (center_dis_t0[v] + center_dis_t1[v])/2

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

        # compere to see if the atom leading the rotation has changed
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
        time steps around when a change in direction event occurs along the x axis

    around_change_co_y : list
        time steps around when a change in direction event occurs along the y axis

    around_change_co_z : list
        time steps around when a change in direction event occurs along the z axis

    times_list : list of arrays
        times when the there is a direction change in x y or z (without time steps around)

    """

    x_vector_array, y_vector_array, z_vector_array = get_step_direction(atoms, md_file, masses)
    vector_arrays = [x_vector_array, y_vector_array, z_vector_array]

    # compare to last step and see if it has changed in x y or z direction
    changes_array_xyz = []
    changes_count_xyz = []
    times_list = []
    for v in range(0, len(vector_arrays)):
        changes_array = []
        for i in range(0, len(x_vector_array) - 1, 1):
            if vector_arrays[v][i] > 0:
                changes_array.append(1)
            elif vector_arrays[v][i] < 0:
                changes_array.append(-1)

    changes_count = 0
    # calculate number of changes
    for j in range(0, len(changes_array) - 1, 1):
        if changes_array[j] != changes_array[j - 1]:
            changes_count = changes_count + 1
            # get times when changes take place
            times_list.append(j)

    changes_count_xyz.append(changes_count)
    changes_array_xyz.append(changes_array)

    # get the positions for the time frame around direction change
    co_xs, co_ys, co_zs, atom_id = trajectories(atoms, md_file)

    around_change_co_x = []
    around_change_co_y = []
    around_change_co_z = []
    for m in range(0, len(co_xs)):
        around_change_x = []
        around_change_y = []
        around_change_z = []
        for x in range(0, len(times_list)):
            # define range around change event
            l_bound = times_list[x] - t_range
            u_bound = times_list[x] + t_range
            # get xyz coordinate for each atom
            for n in range(int(l_bound), int(u_bound)):
                if l_bound >= 0 and u_bound <= len(co_xs[0]):
                    around_change_x.append(co_xs[m][n])
                    around_change_y.append(co_ys[m][n])
                    around_change_z.append(co_zs[m][n])
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
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v):
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    else:
        matrix = np.eye(3)

        return matrix


def expand_points(x_points, y_points, name_points, pos_x, neg_x, pos_y, neg_y, lattice_a, lattice_b, angle):
    """
    Expands points in the x and y direction with respect to hexagonal unit cell

    Parameters
    ----------
    x_points : list
        x coordinates for a set of points

    y_points : list
        y coordinates for a set of points

    name_points : list
        corresponding names of each point

    pos_x : integer
        number of unit cells to move in the positive x direction

    neg_x : integer
        number of unit cells to move in the negative x direction

    pos_y : integer
        number of unit cells to move in the positive y direction

    neg_y : integer
        number of unit cells to move in the negative y direction

    lattice_a : float
        lattice vector b length

    lattice_b : float
        lattice vector b length

    angle : float
        angle between lattice vectors a and b

    Returns
    -------
    new_x_points : list
        x coordinates of added points

    new_y_points : list
        y coordinates of added points

    new_point_names : list
        corresponding names of points

    scale : float
        scale for plotting points with marker="o"
    """
    # calculate project a and b on to x and y (height of cell and additional x correction)
    acute_angel = float(angle) - 90
    y_axis = (math.cos(math.radians(acute_angel))) * lattice_b
    additional_x = -(math.sin(math.radians(acute_angel))) * lattice_b

    new_x_points = []
    new_y_points = []
    new_point_names = []
    for i in range(0, len(x_points)):
        new_x_points_x = []
        new_y_points_x = []
        for j in range(neg_x, pos_x, 1):
            new_x_point = x_points[i] + j * lattice_a
            new_x_points_x.append(new_x_point)
            new_y_points_x.append(y_points[i])
        for k in range(neg_y, pos_y, 1):
            for py in new_y_points_x:
                py = py + k * y_axis
                new_y_points.append(py)
                new_point_names.append(name_points[i])
            for px in new_x_points_x:
                px = px + k * additional_x
                new_x_points.append(px)

    # calculate size of expansion to keep plot points size consistent
    boundary_x = pos_x + abs(neg_x)
    boundary_y = pos_y + abs(neg_y)
    boundary_scale = max(boundary_x, boundary_y)
    scale = 20 / boundary_scale

    # returns 3 arrays of the same length and a scale for plotting size
    return new_x_points, new_y_points, new_point_names, scale


def site_analysis(atoms, md_file, lattice_a, lattice_b, angle, surface_atoms, new_point_dis, new_point_name, time_array):
    """
    Evaluates what surface site an atom is above
    Note some values are hard coded in and may need to be revised depending on system
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
        closest site to the atom

    """

    tx, ty, tz, ta = trajectories(atoms, md_file)
    s_atom_data = get_data(surface_atoms, md_file, 'R')

    s_atom_distances_array = []
    closet_site_array = []
    # for each time step in trajectory find out where the surface atoms are and expand
    for e in range(0, len(time_array), 1):
        t = time_array[e]
        # get atom positions and additional points at time t an put in the same arrays
        # take care with the number of surface atoms to prevent duplication
        adx, ady, names = additional_points(surface_atoms, t, md_file, new_point_dis, new_point_name)
        for k in range(0, len(surface_atoms), 1):
            adx.append(float(s_atom_data[k][t][2]) * 0.529177249)
            ady.append(float(s_atom_data[k][t][3]) * 0.529177249)
            names.append(s_atom_data[k][t][7])
        # find distance from averages
        x_dis_mean_to_point = float(tx[0][t]) - statistics.mean(adx)
        y_dis_mean_to_point = float(ty[0][t]) - statistics.mean(ady)

        acute_angle = float(angle) - 90
        y_axis_height = (math.cos(math.radians(acute_angle))) * lattice_b

        x_dis_ratio = x_dis_mean_to_point/lattice_a
        y_dis_ratio = y_dis_mean_to_point / y_axis_height
        x_dis_ratio_rounded = round(x_dis_ratio)
        y_dis_ratio_rounded = round(y_dis_ratio)

        x_shift_ratio = y_dis_ratio / 1.72  # ratio of x and y shifts may need to change  depending on system
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

        # usefull for checking
        # print(md_file, " time ", t, " min atom distance ", min_s_atom_distance)
        # plt.scatter(x_ps, y_ps, marker="o", s=s)
        # plt.scatter(tx[0][t], ty[0][t], marker="v")
        # plt.axis('square')
        # plt.title(t)
        # plt.show()
        # plt.close()

    return s_atom_distances_array, closet_site_array
