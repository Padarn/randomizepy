import numpy as np
import random
from typing import List, Optional, Union


def check_input_validity(
    N, m, m_unit, m_each, prob, prob_unit, prob_each, num_arms, conditions
):
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer.")
    if m is not None and (not isinstance(m, int) or m < 0 or m > N):
        raise ValueError("m must be a non-negative integer less than or equal to N.")
    if m_unit is not None and (
        not isinstance(m_unit, list)
        or len(m_unit) != N
        or not all(isinstance(x, int) for x in m_unit)
    ):
        raise ValueError("m_unit must be a list of integers with length N.")
    if m_each is not None and (
        not isinstance(m_each, list)
        or len(m_each) != num_arms
        or not all(isinstance(x, int) for x in m_each)
    ):
        raise ValueError("m_each must be a list of integers with length num_arms.")
    if prob is not None and (not isinstance(prob, float) or prob < 0 or prob > 1):
        raise ValueError("prob must be a float between 0 and 1.")
    if prob_unit is not None and (
        not isinstance(prob_unit, list)
        or len(prob_unit) != N
        or not all(isinstance(x, float) for x in prob_unit)
    ):
        raise ValueError("prob_unit must be a list of floats with length N.")
    if prob_each is not None and (
        not isinstance(prob_each, list)
        or len(prob_each) != num_arms
        or not all(isinstance(x, float) for x in prob_each)
    ):
        raise ValueError("prob_each must be a list of floats with length num_arms.")
    if num_arms is not None and (not isinstance(num_arms, int) or num_arms <= 0):
        raise ValueError("num_arms must be a positive integer.")
    if conditions is not None and (
        not isinstance(conditions, list) or len(conditions) != num_arms
    ):
        raise ValueError("conditions must be a list of strings with length num_arms.")


def completely_random_assignment(
    N: int,
    m: Optional[int] = None,
    m_unit: Optional[List[int]] = None,
    m_each: Optional[List[int]] = None,
    prob: Optional[float] = None,
    prob_unit: Optional[List[float]] = None,
    prob_each: Optional[List[float]] = None,
    num_arms: Optional[int] = None,
    conditions: Optional[List[str]] = None,
    check_inputs: bool = True,
) -> List[Union[int, str]]:
    # Basic input checks
    if check_inputs:
        # Add your input validation logic here
        check_input_validity(
            N, m, m_unit, m_each, prob, prob_unit, prob_each, num_arms, conditions
        )

    # Handle the `m_unit` and `prob_unit` cases
    if prob_unit is not None:
        prob = np.unique(prob_unit)[0]
    if m_unit is not None:
        m = np.unique(m_unit)[0]

    # Default conditions for two and multi-arm trials
    if conditions is None:
        conditions = (
            [0, 1]
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


def completely_random_probabilities(
    N: int,
    m: Optional[int] = None,
    m_unit: Optional[List[int]] = None,
    m_each: Optional[List[int]] = None,
    prob: Optional[float] = None,
    prob_unit: Optional[List[float]] = None,
    prob_each: Optional[List[float]] = None,
    num_arms: Optional[int] = None,
    conditions: Optional[List[str]] = None,
    check_inputs: bool = True,
) -> np.ndarray:
    """
    Calculates the probabilities of assignment for a complete random assignment design.

    :param N: The number of units. Must be a positive integer.
    :param m: For a two-arm design, number of units assigned to treatment.
    :param m_unit: For a two-arm design, list of units where each unit should be the same, and length equals N.
    :param m_each: For a multi-arm design, list defining the number of units assigned to each condition.
    :param prob: For a two-arm design, the probability of assigning a unit to treatment.
    :param prob_unit: For a two-arm design, list of probabilities which must be the same for all units.
    :param prob_each: For a multi-arm design, list of probabilities for assigning units to each treatment condition.
    :param num_arms: Number of treatment arms. Inferred from other parameters if not specified.
    :param conditions: Names of the treatment groups.
    :param check_inputs: If True, validates inputs. Defaults to True.
    :return: A NumPy array representing the matrix of probabilities of assignment.
    """

    if check_inputs:
        # Add your input validation logic here
        check_input_validity(
            N, m, m_unit, m_each, prob, prob_unit, prob_each, num_arms, conditions
        )

    # Handle the `m_unit` and `prob_unit` cases
    if prob_unit is not None:
        prob = np.unique(prob_unit)[0]
    if m_unit is not None:
        m = np.unique(m_unit)[0]

    # Default conditions for two and multi-arm trials
    if conditions is None:
        conditions = (
            [0, 1]
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

        prob_matrix = np.zeros((N, 2))
        prob_matrix[:m, 1] = 1
        prob_matrix[m:, 0] = 1
        return prob_matrix

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


# Example usage
Z = completely_random_assignment(N=100, m=50)
print(np.bincount(Z))

# This is a direct translation. Some Python-specific optimizations or stylistic changes might be needed.
