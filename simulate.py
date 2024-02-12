import ray
from completely_random_assignment import completely_random_assignment

@ray.remote
def remote_completely_random_assignment(N, m):
    return completely_random_assignment(N, m)


def simulate(n_sim=100):
    futures = [remote_completely_random_assignment.remote(100, 50)]
    results = [ray.get(f) for f in futures]
    return results

if __name__ == '__main__':
    ray.init()
    sim = simulate()
    print(sim)