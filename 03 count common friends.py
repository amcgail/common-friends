from common import *

def sample_until_common_friend(net, k=2):
    """
    Sample random individuals until we find a "common friend to k" --
    someone who appears in the friend lists of k different sampled individuals.
    
    This implements the sampling procedure from "The Power of Common Friends":
    we draw random people and track how many sampled individuals list each
    person as a friend. We stop when some person has been listed k times.
    
    Args:
        net: Network object with .nodes list and .nbds dict (node -> set of neighbors)
        k: Number of times someone must appear in friend lists (default 2)
    
    Returns:
        (sample_size, common_friend_id): Number of people sampled, and the ID of 
        a randomly chosen common friend to k (if multiple exist with same count)
    """
    # Count how many sampled individuals list each person as a friend
    friend_counts = defaultdict(int)
    
    # Track which nodes we've already sampled
    sampled = set()
    
    # Buffer of candidates for efficient random sampling
    buffer = None
    
    while not friend_counts or max(friend_counts.values()) < k:
        # Refill buffer if empty
        if buffer is None or len(buffer) == 0:
            candidates = [n for n in net.nodes if n not in sampled]
            buffer = list(np.random.choice(candidates, min(200, len(candidates)), replace=False))
        
        # Sample a random individual
        sampled_node = buffer.pop()
        sampled.add(sampled_node)
        
        # Update friend counts for all neighbors of the sampled individual
        for neighbor in net.nbds[sampled_node]:
            friend_counts[neighbor] += 1
    
    # Find all individuals who achieved the max count (common friends to k)
    max_count = max(friend_counts.values())
    candidates = [node for node, count in friend_counts.items() if count == max_count]
    common_friend = np.random.choice(candidates)
    
    return len(sampled), common_friend


def run_simulation_study(net, ks=(1, 2, 3, 4, 5, 6), n_trials=300):
    """
    Run the common friends simulation study for a network.
    
    For each value of k, we repeatedly sample until finding a common friend to k,
    and record both the sample size required and the degree of the common friend found.
    
    Args:
        net: Network object
        ks: Tuple of k values to test (number of common friends required)
        n_trials: Number of simulation trials per k value
    """
    print(f'Starting simulation for {net.name}')
    
    # Precompute degree statistics for the network
    degrees = net.ds
    max_degree_in_network = degrees.max()
    top_5_degrees = set(sorted(degrees)[-5:])
    
    results = defaultdict(list)
    
    for k in ks:
        print(f'  k={k}: sampling until common friend to {k}...')
        
        found_degrees = []      # Degree of the common friend found in each trial
        sample_sizes = []       # Number of individuals sampled in each trial
        
        for _ in range(n_trials):
            sample_size, common_friend = sample_until_common_friend(net, k=k)
            degree = net.ddict[common_friend]
            
            found_degrees.append(degree)
            sample_sizes.append(sample_size)
        
        found_degrees = np.array(found_degrees)
        sample_sizes = np.array(sample_sizes)
        
        # Record statistics
        results["avg_deg"].append(np.mean(found_degrees))
        results["med_deg"].append(np.median(found_degrees))
        results["min_deg"].append(np.min(found_degrees))
        results["max_deg"].append(np.max(found_degrees))
        results["5pct_deg"].append(np.percentile(found_degrees, 5))
        results["95pct_deg"].append(np.percentile(found_degrees, 95))
        
        results["avg_samp"].append(np.mean(sample_sizes))
        results["5pct_samp"].append(np.percentile(sample_sizes, 5))
        results["95pct_samp"].append(np.percentile(sample_sizes, 95))
        
        results["got_max"].append(np.mean(found_degrees == max_degree_in_network))
        results["got_top5"].append(np.mean([d in top_5_degrees for d in found_degrees]))
    
    # Write results to file
    with open(f'tables/{net.name}.sim.txt', 'w') as f:
        for i, k in enumerate(ks):
            f.write(f'------ k={k} ------\n')
            for stat_name, stat_values in results.items():
                f.write(f'{stat_name}: {stat_values[i]}\n')
            
            # What percentile is the average degree found?
            avg_deg = results["avg_deg"][i]
            percentile = 1 - (degrees >= avg_deg).sum() / len(degrees)
            f.write(f'avg_deg percentile: {percentile:0.3%}\n')


if __name__ == '__main__':
    run_simulation_study(NOLA)
    run_simulation_study(citations)