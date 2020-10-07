# import math
# import random
# import statistics


# def euclidean_distance_v(vector1, vector2):
#     '''
#     Calculate the euclidean distance between
#     two vectors.
#     '''
#     value = 0
#     for v1, v2 in zip(vector1, vector2):
#         value += (v1 - v2)**2
#     return math.sqrt(value)


# def multiple_euclidean_distance_v(vlist1, vlist2):
#     '''
#     Extending the above function
#     '''
#     dist_list = []
#     for list1, list2 in zip(vlist1, vlist2):
#         value = 0
#         for v1, v2 in zip(list1, list2):
#             value += (v1 - v2)**2
#         dist_list.append(value)
#     mean_dist = sum(dist_list) / len(dist_list)
#     var_dist = statistics.variance(dist_list)
#     return mean_dist, var_dist


# # Function to generate a list of lists containing frequencies between 0 and 1
# # and summing up to 1
# def rand_frequencies(n, length):
#     '''
#     List of lists containing random frequencies
#     N is how many lists and length is length of each list.
#     '''
#     return_list = []
#     for i in range(n + 1):
#         current_list = []
#         for j in range(length + 1):
#             current_list.append(random.uniform(0, 1))
#         append_values = [(val / sum(current_list)) for val in current_list]
#         return_list.append(append_values)
#     return return_list


# a = rand_frequencies(10, 5)
# b = rand_frequencies(10, 5)

# data = multiple_euclidean_distance_v(a, b)
# print(data)


import subprocess
import os

command = f"esearch -db assembly -query GCA_900066305.1 | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
get_download_link = subprocess.Popen(
    command, stdout=subprocess.PIPE,
    shell=True)
ftp_link = get_download_link.stdout.read().strip().decode("utf-8")
file_link = f"{str(os.path.basename(ftp_link))}_genomic.fna.gz"
full_link = f"{ftp_link}/{file_link}"
print(f"FTP link: {ftp_link}")
print(f"File link: {file_link}")
print(f"Full link: {full_link}")
