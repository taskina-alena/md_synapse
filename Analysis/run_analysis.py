from Oct2023.Analysis.cluster_analysis import ClusterAnalysis
import os
import time

if __name__ == '__main__':

	eps = 6
	num_vesicles = 200
	syn_density = 0.2
	model_type = 'mono'
	if model_type=='bi':
		bi_domain = True
	else:
		bi_domain = False

	file_pattern = f'2d_{model_type}_eps{eps}_nv{num_vesicles}_sd{syn_density}'
	scratch_loc = 'Results/'

	movie_name = scratch_loc + 'Movie_%s.xyz' % file_pattern
	oligo_name = scratch_loc + 'Oligo_%s.dat' % file_pattern

	print(movie_name)
	print(os.path.exists(movie_name))
	frames = 'all'
	start_time = time.time()
	if os.path.exists(movie_name):
		va = ClusterAnalysis(movie_name, num_vesicles, bi_domain)
		va.clusters_over_time(oligo_name, frames)

	print("--- %s seconds ---" % (time.time() - start_time))

	va.plot_clusters(-1)


