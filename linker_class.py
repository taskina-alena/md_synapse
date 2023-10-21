import numpy as np

class LinkedVesicles:
    """Generate linked vesicles with specified attributes."""

    def __init__(self, num_v, nums_syn, side, radii, bi_domain, kBend, kBond_vesice, kBond_itra_syn, eps):
        """
        Initialize linked vesicles object.

        :param num_v: Number of vesicles.
        :param nums_syn: Number of synapsins for each vesicle.
        :param side: Size of the box side.
        :param radii: Dictionary containing radii for vesicle and synapsin.
        :param bi_domain: Boolean indicating if bi-domain structure is used.
        """
        # Initialization of basic attributes
        self.bi_domain = bi_domain
        self.nums_syn = nums_syn
        self.num_v = int(num_v)
        self.side = side
        self.radii_vesicles = radii['vesicles']
        self.radius_syn = radii['synapsin']
        self.lj_factor = 2 ** (1 / 6)
        self.kBend = kBend
        self.kBond_vesicle = kBond_vesice
        self.kBond_intra_syn = kBond_itra_syn
        self.eps = eps


        # Initialize internal data structures and parameters
        self._initialize_data_structures()

        # Initialize lattice structure and setup vesicles
        self._initialize_lattice()
        self.build_vesicles()

        if self.bi_domain:
            self.num_types += 1

        self._setup_types_and_masses()

        if self.bi_domain:
            self._describe_angle()

        self._describe_bonds()
        self._describe_interactions()

    def _initialize_data_structures(self):
        """Initialize data structures for coordinates, bonds, and angles."""
        # Lists to hold definitions of atoms, bonds, and angles respectively
        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.bonds_interactions = ["bond_style harmonic\n"]
        self.interactions = [f"pair_style cosine/squared 3 \n"]
        # Initialize counters
        self.numAll, self.c, self.num_bonds, self.bond_id, self.num_angles, self.angle_id, self.num_types = (0,) * 7
        # List to hold types and their associated masses
        self.type_mass_list = []
        if self.bi_domain:
            self.angles = ["\nAngles \n \n"]
            self.angles_interactions = ['angle_style harmonic\n']

    def _initialize_lattice(self):
        """Setup the 2D lattice for vesicle positioning."""
        # Dimensions of the system box
        self.Lx, self.Ly, self.Lz = self.side, self.side, 0.1
        # Compute lattice constant based on Lennard-Jones potential factor
        self.lattice_constant = 2 * (self.lj_factor * max(self.radii_vesicles) + 2 * self.radius_syn)
        # Compute lattice sites in 2D
        self.lattice_sites = self._init_lattice_2d()
        # Indices for the lattice sites
        self.lattice_inds = np.arange(self.lattice_sites.shape[0])
        # Randomly select starting indices for vesicles
        self.start_inds_v = np.random.choice(self.lattice_inds, size=self.num_v, replace=False)
        self.ves_lattice_sites = self.lattice_sites[self.start_inds_v]

    def _init_lattice_2d(self):
        """Initialize a 2D lattice of points."""
        # Number of lattice points in x and y direction
        num_x, num_y = int(self.Lx / self.lattice_constant), int(self.Ly / self.lattice_constant)
        # Create an array of x and y coordinates
        x_coords = np.linspace(-num_x * self.lattice_constant / 2, num_x * self.lattice_constant / 2, num_x)
        y_coords = np.linspace(-num_y * self.lattice_constant / 2, num_y * self.lattice_constant / 2, num_y)
        # Create a mesh grid for x and y coordinates
        xx, yy = np.meshgrid(x_coords[1:], y_coords[1:])
        # Return coordinates in a 2D array format
        return np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

    def _dist_on_circle(self, radius, num_syn):
        """Calculate positions on the circumference of a circle for placing synapsins."""
        # Calculate circumference of the circle
        cf = 2 * np.pi * radius
        # Check if the number of synapsins is feasible for the given radius
        if num_syn > np.floor(cf / self.radius_syn):
            raise RuntimeError(f'Too many linkers on circumference for radius {radius}. Reduce synapsin density or adjust sigma.')
        # Compute the angle between each synapsin on the circle
        if num_syn == 0:
            return []
        else:
            linker_angle = 2.0 * np.pi / num_syn
            angles = np.arange(num_syn) * linker_angle
            # Return x and y coordinates of synapsins
            return np.column_stack([radius * np.sin(angles), radius * np.cos(angles)])

    def _setup_types_and_masses(self):
        """Initialize types and their masses for the system."""
        self.type_mass_list = [[i + 1, 1] for i in range(self.num_types)]

    def build_vesicles(self):
        """Construct vesicles based on provided parameters."""
        # For each vesicle, add its center and its surrounding synapsins
        for ves, (x, y, z) in enumerate(self.ves_lattice_sites):
            self.c += 1
            radius, num_syn = self.radii_vesicles[ves], self.nums_syn[ves]
            self._add_central_vesicle(radius, ves + 1, x, y, z)
            inner_radius = (radius + self.radius_syn) * self.lj_factor
            outer_radius = inner_radius + 2 * self.radius_syn
            if self.bi_domain:
                self._add_shell_syn(num_syn, inner_radius, x, y, z, False, self.numAll)
                self._add_shell_syn(num_syn, outer_radius, x, y, z, True, self.numAll)
            else:
                self._add_shell_syn(num_syn, radius * self.lj_factor, x, y, z, True, self.numAll)

    def _add_central_vesicle(self, radius, particle_type, x, y, z):
        """Add the central vesicle particle to the coordinates list."""
        self.numAll += 1
        self.num_types += 1
        coord_string = f"\t {self.numAll} {self.c} {particle_type} {radius} {x} {y} {z} 0 0 0 \n"
        self.coords.append(coord_string)

    def _add_shell_syn(self, num_syn, radius, x, y, z, outer, ves_id):
        """Add synapsins around the vesicle, either in the inner shell or outer shell, outer shell lacks intra interactions"""
        particle_type = self.num_v + (1 if self.bi_domain else 0) + (self.c if outer else 0)
        syn_positions = self._dist_on_circle(radius, num_syn)
        # Add each synapsin to the coordinates list
        for pos in syn_positions:
            self._add_syn_to_shell(pos, x, y, z, particle_type, num_syn, outer, ves_id)
        if outer:
            self.num_types += 1

    def _add_syn_to_shell(self, pos, x, y, z, particle_type, num_syn, outer, ves_id):
        """Add individual synapsin to the shell."""

        # add coordinates
        x_shell, y_shell = x + pos[0], y + pos[1]
        self.numAll += 1
        coord_string = f"\t {self.numAll} {self.c} {particle_type} {self.radius_syn} {x_shell} {y_shell} {z} 0 0 0 \n"
        self.coords.append(coord_string)

        # add bond
        self.bond_id += 1
        self.num_bonds += 1
        if (outer and self.bi_domain):
            bond_target = self.numAll - num_syn
            bond_type = self.num_v +1
        else:
            bond_target = ves_id
            bond_type = self.c
        bond_string = f"\t {self.bond_id} {bond_type} {self.numAll} {bond_target} \n"
        self.bonds.append(bond_string)

        # add angle
        if outer and self.bi_domain:
            self.angle_id += 1
            self.num_angles += 1
            angle_string = f"\t {self.angle_id} 1 {self.c} {self.numAll - num_syn} {self.numAll}\n"
            self.angles.append(angle_string)

    def _describe_bonds(self):
        if self.bi_domain:
            self.bonds_interactions.append(f"bond_coeff {self.num_v+1:n} {self.kBond_intra_syn} {2*self.radius_syn}\n")
        for i in range(self.num_v):
            radius = self.lj_factor * (self.radii_vesicles[i] + self.radius_syn)
            self.bonds_interactions.append(f"bond_coeff {i+1:n} {self.kBond_vesicle} {radius}\n")

    def _describe_angle(self):
        ANGLE = 180
        self.angles_interactions.append(f"angle_coeff 1 {self.kBend} {ANGLE}\n")

    def _add_interaction(self, type1, type2, r1, r2, attr=False):
        cutoff = self.lj_factor*(r1+r2)
        if attr:
            cutoff += 1
        self.interactions.append(f"pair_coeff {type1:n} {type2:n} {self.eps} {self.lj_factor*(r1+r2)} {cutoff} wca \n")

    def _describe_interactions(self):
        #repulsive interactions of vesicles with other vesicles and synapsin
        for i in range(self.num_v):
            i_radius = self.radii_vesicles[i]
            for j in range(self.num_v):
                j_radius = self.radii_vesicles[j]
                if i <= j:
                    self._add_interaction(i+1, j+1, i_radius, j_radius)
                self._add_interaction(i+1, self.num_v+j+1, i_radius, self.radius_syn)
            if self.bi_domain:
                self._add_interaction(i+1, 2*self.num_v+1, i_radius, self.radius_syn)

        #repulsive interaction of inner shell within itself and with other synapsins
        if self.bi_domain:
            self._add_interaction(self.num_v+1, self.num_v+1, self.radius_syn, self.radius_syn)
            for j in range(self.num_v):
                self._add_interaction(self.num_v+1, self.num_v+1+j+1, self.radius_syn, self.radius_syn)

        #interactions within outer shell
        offset = self.num_v + 1 if self.bi_domain else self.num_v

        for i in range(self.num_v):
            for j in range(self.num_v):
                if i<j:
                    self._add_interaction(offset+i+1, offset+j+1, self.radius_syn, self.radius_syn, attr=True)
                elif i==j:
                    self._add_interaction(offset + i + 1, offset + j + 1, self.radius_syn, self.radius_syn, attr=False)






