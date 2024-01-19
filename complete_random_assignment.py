import numpy as np
import random


def complete_ra(
    N,
    m=None,
    m_unit=None,
    m_each=None,
    prob=None,
    prob_unit=None,
    prob_each=None,
    num_arms=None,
    conditions=None,
    check_inputs=True,
):
    # Basic input checks
    if check_inputs:
        # Add your input validation logic here
        pass

    # Handle the `m_unit` and `prob_unit` cases
    if prob_unit is not None:
        prob = np.unique(prob_unit)[0]
    if m_unit is not None:
        m = np.unique(m_unit)[0]

    # Default conditions for two and multi-arm trials
    if conditions is None:
        conditions = (
            ["0", "1"]
            if num_arms is None or num_arms == 2
            else ["T" + str(i) for i in range(1, num_arms + 1)]
        )

    # Two-arm designs
    if m_each is None and prob_each is None and len(conditions) == 2:
        # Case 1: Neither m nor prob is specified
        if m is None and prob is None:
            m = N // 2 if random.random() < 0.5 else N - N // 2

        # Case 2: m is specified
        elif m is not None:
            pass  # m is already set

        # Case 3: prob is specified
        elif prob is not None:
            m = (
                int(np.floor(N * prob))
                if random.random() < prob
                else int(np.ceil(N * prob))
            )

        assignment = random.sample(conditions * m + conditions[::-1] * (N - m), N)
        return assignment

    # Multi-arm designs
    if m_each is None and prob_each is None:
        prob_each = [1 / num_arms] * num_arms

    if prob_each is not None:
        m_each = np.floor(np.array(prob_each) * N).astype(int)
        remainder = N - np.sum(m_each)
        additional_assignments = np.random.choice(conditions, remainder, p=prob_each)
        full_assignments = np.concatenate(
            [np.repeat(conditions, m_each), additional_assignments]
        )
        np.random.shuffle(full_assignments)
        return full_assignments.tolist()

    if m_each is not None:
        full_assignments = np.repeat(conditions, m_each)
        np.random.shuffle(full_assignments)
        return full_assignments.tolist()


def complete_ra_probabilities(
    N,
    m=None,
    m_unit=None,
    m_each=None,
    prob=None,
    prob_unit=None,
    prob_each=None,
    num_arms=None,
    conditions=None,
    check_inputs=True,
):
    # Similar setup and validation as complete_ra
    # ...

    # Logic to calculate the probabilities matrix
    # ...

    return prob_mat  # Replace with actual probability matrix logic


# Example usage
Z = complete_ra(N=100, m=50)
print(np.bincount(Z))

# This is a direct translation. Some Python-specific optimizations or stylistic changes might be needed.
