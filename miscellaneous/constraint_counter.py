GATES_PER_HASH = 2100
RHO = 16
COLLAPSING_FACTOR = 8
LAST_LEVEL_SIZE = 1
DOMAIN_SIZE = 2^24


def merkle_proof_num_cnstr(tree_size):
    # tree size is actually the same as num of leaves
    # we are using Blake2s - so there is no padding, and 
    # 2:1 transform will take GATES_PER_HASH
    # there will be log(tree_size) hashes in the Merle path
    num_hashes = log(tree_size, 2)
    return GATES_PER_HASH * num_hashes


def calculate_num_cnstrs_per_query():
    total_cnsts_count = 0
    # all elements of the same coset are located in the single leaf of our merklee tree 
    # so the total number of leaves on the top level will be DOMAIN_SIZE / COLLAPSING_FACTOR
    # and with each level this number will be reduced by COLLAPSING_FACTOR
    cur_level_size = DOMAIN_SIZE * RHO
    while cur_level_size > LAST_LEVEL_SIZE:
        num_of_leaves = cur_level_size / COLLAPSING_FACTOR
        total_cnsts_count += merkle_proof_num_cnstr(num_of_leaves)
        cur_level_size /= COLLAPSING_FACTOR
    return total_cnsts_count


NUM_OF_QUERIES = 16

def calculate_num_cnsts_per_cicruit():
    # wire polynomials a, b, c, d are packed in the single tree
    # t_0, t_mid, t_high are packed in the single tree
    # z_1, z_2 are packed in the same tree
    # all selectors are also packed in the single tree
    # so, there will be 4 trees in total
    return 4 * NUM_OF_QUERIES * calculate_num_cnstrs_per_query()


print calculate_num_cnstrs_per_query()
print calculate_num_cnsts_per_cicruit()