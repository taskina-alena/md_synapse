from tools import lmp_tools as lt
import freud as fr
import numpy as np
from matplotlib import pyplot as plt
import time
start_time = time.time()
from scipy.spatial import ConvexHull

class ClusterAnalysis(object):

    def __init__(self, movie_name, num_vesicles, bi_domain):
        self.movie_name = movie_name
        self.num_vesicles = num_vesicles
        self.bi_domain = bi_domain

        self.r_cut = 2.5 * 1  ## cut-off for cluster computes
        movie_format = ['id', 'type', 'mol', 'x', 'y', 'z', 'radius']  # structure of the movie data to be loaded

        self.max_k_mer = self.num_vesicles
        movie_format = ['id', 'type', 'mol', 'x', 'y', 'z', 'radius']  # structure of the movie data to be loaded

        read_lines = lt.header_lines(self.movie_name)  ## select or slice from read_lines to load single frames
        self.work_movie = lt.split_movie(self.movie_name, movie_format, read_lines)

        box_data = np.genfromtxt(self.movie_name, skip_header=5, max_rows=3)  # grab box data from movie file
        box_lens = box_data[:, 1] - box_data[:, 0]
        self.box = fr.box.Box(box_lens[0], box_lens[1], is2D=True)  # define the box that the freud package can use

    def _parse_config(self,):
        all_types = np.sort(np.unique(self.config[:, 1]))  # sorted list of all types appearing in the simualtion
        types_vesicles = all_types[:self.num_vesicles]
        types_synapsin = all_types[-self.num_vesicles:]  # Number of types can vary (e.g. if vesicles are unoccupied)
        self.lines_vesicles = np.where(np.isin(self.config[:, 1], types_vesicles))[0]
        self.lines_synapsin = np.where(np.isin(self.config[:, 1], types_synapsin))[0]

    def _compute_syn_clusters(self):
        # Setup the neighbour search based on locations of the synapsin tails
        coords_synapsin = self.config[self.lines_synapsin][:,3:6]
        aq = fr.AABBQuery(self.box, coords_synapsin)
        query_result = aq.query(coords_synapsin, dict(r_max=self.r_cut))
        n_list = query_result.toNeighborList()
        return n_list

    def _shift_from_syn_to_ves(self, n_list):
        # Transfer the neighbour list from index to types
        config_synapsin = self.config[self.lines_synapsin]
        n_list_syn_types = np.transpose(np.vstack([config_synapsin[n_list[:, 0], 1], config_synapsin[n_list[:, 1], 1]]))
        n_list_syn_types_unique = n_list_syn_types[n_list_syn_types[:, 0] != n_list_syn_types[:, 1]]
        if self.bi_domain:
            n_list_ves = n_list_syn_types_unique - 1 - self.num_vesicles
        else:
            n_list_ves = n_list_syn_types_unique - self.num_vesicles
        return np.unique(n_list_ves, axis=0)

    def _ves_clusters_props(self, n_list):

        config_ves = self.config[self.lines_vesicles]
        config_ves_type_sorted = config_ves[config_ves[:, 1].argsort()]  # Permute the config such that it is sorted by types index
        coords_ves_type_sorted = config_ves_type_sorted[:, 3:6]

        # Transfer info into neighbour list object
        query_point_indices = n_list[:, 0].astype(int) - 1
        point_indices = n_list[:, 1].astype(int) - 1
        distances = self.box.compute_distances(coords_ves_type_sorted[query_point_indices],coords_ves_type_sorted[point_indices])
        num_query_points = len(coords_ves_type_sorted)
        num_points = len(coords_ves_type_sorted)
        nlist = fr.locality.NeighborList.from_arrays(num_query_points, num_points, query_point_indices, point_indices, distances)

        #create cluster from neighbor list
        aq = fr.AABBQuery(self.box, coords_ves_type_sorted)
        cl = fr.cluster.Cluster()
        cl.compute(aq, neighbors=nlist)
        cl_prop = fr.cluster.ClusterProperties()
        cl_prop.compute(aq, cl.cluster_idx)

        # Extract areas and perimeters for each cluster
        cluster_areas = []
        cluster_perimeters = []
        for cluster_id in range(cl.num_clusters):
            cluster_points = coords_ves_type_sorted[cl.cluster_keys[cluster_id]]
            area, perimeter = self._compute_cluster_properties(cluster_points[:, :2])  # Only x, y coordinates
            cluster_areas.append(area)
            cluster_perimeters.append(perimeter)

        return cl, cl_prop, cluster_areas, cluster_perimeters

    def _find_ves_clusters(self, frame):

        self._parse_config()
        n_list_syn = self._compute_syn_clusters()
        n_list_ves = self._shift_from_syn_to_ves(n_list_syn)
        return self._ves_clusters_props(n_list_ves)


    def _save_cl_props(self, frame, cl_prop, cluster_areas, cluster_perimeters):
        mean_oligo_size = np.mean(cl_prop.sizes)

        # Calculate circularity for each cluster and then compute the mean
        circularity = [4 * np.pi * area / (perimeter ** 2) for area, perimeter in
                       zip(cluster_areas, cluster_perimeters) if area]
        mean_circularity = np.mean(circularity)

        oligo_dist = np.zeros(shape=self.max_k_mer)
        size_dist = np.unique(cl_prop.sizes, return_counts=True)
        np.testing.assert_array_equal(size_dist[0], np.sort(size_dist[0]))
        assert np.max(size_dist[0]) <= self.max_k_mer
        for size, count in zip(size_dist[0], size_dist[1]):
            oligo_dist[size] = count

        self.oligo_dist_time[frame, 0] = frame
        self.oligo_dist_time[frame, 1] = mean_circularity
        self.oligo_dist_time[frame, 2] = mean_oligo_size
        self.oligo_dist_time[frame, 3:] = oligo_dist

    def clusters_over_time(self, oligo_file_name, frames_to_process):

        if frames_to_process=='all':
            frames_to_process = np.arange(len(self.work_movie))

        self.oligo_dist_time = np.zeros(shape=(len(frames_to_process), self.max_k_mer + 3))

        for frame in frames_to_process:  # main loop through simulation time frames

            self.config = self.work_movie[frame] # grab the full configuration at time 'frame'
            cl, cl_prop, areas, perimeters = self._find_ves_clusters(frame)
            self._save_cl_props(frame, cl_prop, areas, perimeters)

        np.savetxt(oligo_file_name, self.oligo_dist_time, delimiter=' ')

    def _compute_cluster_properties(self, cluster_points):
        """
        Compute the area and perimeter of the cluster's convex hull.
        """
        if len(cluster_points)>2:
            hull = ConvexHull(cluster_points)
            area = hull.volume
            perimeter = hull.area
        else:
            area = None
            perimeter = None

        return area, perimeter

    def plot_clusters(self, frame):

        self.config = self.work_movie[frame]
        cl, cl_prop, area, perimeter = self._find_ves_clusters(frame)

        # Extract vesicle coordinates for plotting
        config_ves = self.config[self.lines_vesicles]
        config_ves_type_sorted = config_ves[config_ves[:, 1].argsort()]  # Permute the config such that it is sorted by types index
        coords_ves = config_ves_type_sorted[:, 3:6]

        # Create a new figure and axis for plotting
        fig, ax = plt.subplots(figsize=(9, 6))

        # Color each cluster with a different color
        cluster_colors = plt.cm.jet(np.linspace(0, 1, cl.num_clusters))

        # Plot each cluster with its unique color
        for cluster_id, color in enumerate(cluster_colors):
            cluster_points = coords_ves[cl.cluster_keys[cluster_id]]
            ax.scatter(cluster_points[:, 0], cluster_points[:, 1], color=color, label=f"Cluster {cluster_id}")

        # Set plot title, labels, and legend
        ax.set_title(f"Vesicle Clusters for Frame {frame}, " + 'mean cluster size ' + str(np.mean(cl_prop.sizes)))
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")

        # Display the plot
        plt.tight_layout()
        plt.show()




