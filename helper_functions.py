def check_randomizr_arguments(
    N=None,
    prob=None,
    prob_unit=None,
    m=None,
    m_unit=None,
    m_each=None,
    prob_each=None,
    blocks=None,
    block_m=None,
    block_m_each=None,
    block_prob=None,
    block_prob_each=None,
    clusters=None,
    num_arms=None,
    simple=None,
    conditions=None,
    **kwargs,
):
    if clusters is not None:
        N = len(set(clusters))

    if blocks is not None:
        if simple:
            raise ValueError(
                "You can't specify `simple` when using blocked assignment."
            )

        if clusters is not None:
            N_per_block = {b: len(set(c)) for b, c in zip(blocks, clusters)}
            if any(sum(table(blocks, clusters) > 0) > 1):
                raise ValueError(
                    "All units within a cluster must be in the same block."
                )
        else:
            N_per_block = {b: blocks.count(b) for b in blocks}

        if N is None:
            N = sum(N_per_block.values())
        elif N != sum(N_per_block.values()):
            raise ValueError("N should equal the length of blocks.")

        N_blocks = len(N_per_block)
        blocks = {"N_per_block": N_per_block, "N_blocks": N_blocks}

    if N is None:
        raise ValueError("N, blocks or clusters must be specified.")

    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer scalar.")

    if len(set(conditions)) != len(conditions):
        raise ValueError("You must supply unique values to conditions.")

    if num_arms is not None and conditions is not None and num_arms != len(conditions):
        raise ValueError(
            "If both conditions and num_arms are specified, the length of conditions must be equal to num_arms."
        )

    conflict_args = [
        "prob",
        "prob_each",
        "prob_unit",
        "m",
        "m_unit",
        "m_each",
        "block_prob",
        "block_prob_each",
        "block_m",
        "block_m_each",
    ]
    specified_args = {
        arg: value for arg, value in kwargs.items() if arg in conflict_args
    }

    if len(specified_args) > 1:
        raise ValueError(
            "Please specify only one of " + " and ".join(specified_args.keys()) + "."
        )
    elif len(specified_args) == 1:
        arg = next(iter(specified_args))
        arg_block = arg.startswith("block_")
        arg_each = arg.endswith("_each")

        if arg_block and blocks is None:
            raise ValueError("Specified `" + arg + "` but blocks is None.")

        if simple and "prob" not in arg:
            raise ValueError("You can't specify `" + arg + "` when simple = True.")

        if num_arms is not None and conditions is not None:
            _check_ra_arg_num_arms_conditions(
                arg, arg_block, arg_each, specified_args[arg], num_arms, conditions
            )

        _check_ra[arg](
            N, blocks, clusters, num_arms, conditions, simple, specified_args[arg]
        )

    if num_arms is None:
        if conditions is not None:
            num_arms = len(conditions)
        elif len(specified_args) == 0:
            num_arms = 2
        elif not arg_each:
            num_arms = 2
        elif not arg_block and not isinstance(specified_args[arg], list):
            num_arms = len(specified_args[arg])
        else:
            num_arms = len(specified_args[arg])

        if num_arms == 2 and conditions is None:
            conditions = [0, 1]

    if conditions is None:
        conditions = ["T" + str(i) for i in range(1, num_arms + 1)]

    ret = {
        "num_arms": num_arms,
        "conditions": conditions,
        "N_per_block": blocks.get("N_per_block"),
    }

    return ret


class CheckRA:
    def __init__(self):
        conflict_args = [
            "prob",
            "prob_each",
            "prob_unit",
            "m",
            "m_unit",
            "m_each",
            "block_prob",
            "block_prob_each",
            "block_m",
            "block_m_each",
        ]

    @classmethod
    def check_arguments(cls, **kwargs):
        specified_args = {
            k: v for k, v in kwargs.items() if k in cls.conflict_args and v is not None
        }

        if len(specified_args) > 1:
            raise ValueError(
                f"Please specify only one of {', '.join(specified_args.keys())}."
            )
        elif len(specified_args) == 1:
            arg, value = next(iter(specified_args.items()))
            arg_block = arg.startswith("block_")
            arg_each = arg.endswith("_each")

            if arg_block and kwargs.get("blocks") is None:
                raise ValueError(f"Specified `{arg}` but blocks is None.")

            if kwargs.get("simple") is True and "prob" not in arg:
                raise ValueError(f"You can't specify `{arg}` when simple = True.")

            getattr(cls, f"check_{arg}")(
                arg,
                arg_block,
                arg_each,
                value,
                kwargs.get("num_arms"),
                kwargs.get("conditions"),
            )
            getattr(cls, arg)(
                kwargs.get("N"),
                kwargs.get("blocks"),
                kwargs.get("clusters"),
                kwargs.get("num_arms"),
                kwargs.get("conditions"),
                kwargs.get("simple"),
                value,
            )

    def check_ra_arg_num_arms_conditions(
        arg, arg_block, arg_each, value, num_arms, conditions
    ):
        w = 2 if not arg_each else (len(value) if not arg_block else value.shape[1])

        if num_arms is not None and num_arms != w:
            num_arms_suffix = (
                " and num_arms must be 2"
                if not arg_each
                else (" equal to num_arms" if not arg_block else " equal to num_arms")
            )
            raise ValueError(
                f"If {arg} and num_arms are both specified, {arg} must be length {w}{num_arms_suffix}."
            )

        if conditions is not None and len(conditions) != w:
            conditions_suffix = (
                "" if not arg_each else " equal to the length of conditions"
            )
            raise ValueError(
                f"If {arg} and conditions are both specified, {arg} and conditions must both be length {w}{conditions_suffix}."
            )

    def prob(N, blocks, clusters, num_arms, conditions, simple, prob):
        if any(p > 1 or p < 0 for p in prob):
            raise ValueError(
                "The probability of assignment to treatment must be between 0 and 1."
            )
        if len(prob) not in [1]:
            raise ValueError("`prob` must be of length 1.")

    def prob_unit(self, N, blocks, clusters, num_arms, conditions, simple, prob_unit):
        if any(p > 1 or p < 0 for p in prob_unit):
            raise ValueError(
                "The probability of assignment to treatment must be between 0 and 1."
            )

        if blocks is None:
            if not simple:
                if not all(p == prob_unit[0] for p in prob_unit):
                    raise ValueError(
                        "In a complete random assignment design, `prob_unit` must be the same for all units."
                    )
        else:
            if not all(p == prob_unit[0] for p in prob_unit):
                raise ValueError(
                    "In a block random assignment design, `prob_unit` must be the same for all units within the same block."
                )

        if clusters is not None:
            if len(prob_unit) != len(clusters):
                raise ValueError("`prob_unit` must be of length N.")
        else:
            if len(prob_unit) != N:
                raise ValueError("`prob_unit` must be of length N.")

    def check_ra_prob_each(
        self, N, blocks, clusters, num_arms, conditions, simple, prob_each
    ):
        if any(p > 1 or p < 0 for p in prob_each):
            raise ValueError(
                "The probabilities of assignment to any condition may not be greater than 1 or less than zero."
            )

        if isinstance(prob_each, list):
            if sum(prob_each) != 1:
                raise ValueError(
                    "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs."
                )
        elif isinstance(prob_each, np.ndarray):
            if any(np.sum(prob_each, axis=1) != 1):
                raise ValueError(
                    "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs."
                )
            if prob_each.shape[0] not in [1, N]:
                raise ValueError("`prob_each` must have either 1 or N rows.")
        else:
            raise ValueError("`prob_each` must be a list or numpy array.")

    def m(self, N, blocks, clusters, num_arms, conditions, simple, m):
        if any(p > 1 or p < 0 for p in prob_each):
            raise ValueError(
                "The probabilities of assignment to any condition may not be greater than 1 or less than zero."
            )
        if isinstance(prob_each, list):
            if sum(prob_each) != 1:
                raise ValueError(
                    "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs."
                )
        elif isinstance(prob_each, np.ndarray):
            if any(np.sum(prob_each, axis=1) != 1):
                raise ValueError(
                    "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs."
                )

    def m_unit(self, N, blocks, clusters, num_arms, conditions, simple, m_unit):
        # if it's complete
        if blocks is None:
            if not all(x == m_unit[0] for x in m_unit):
                raise ValueError(
                    "In a complete random assignment design, `m_unit` must be the same for all units."
                )
        else:
            # if it's blocked
            for block in blocks:
                if not all(x == m_unit[block[0]] for x in m_unit[block]):
                    raise ValueError(
                        "In a block random assignment design, `m_unit` must be the same for all units within the same block."
                    )

        # if it's clustered:
        if clusters is not None:
            if len(m_unit) != len(clusters):
                raise ValueError("`m_unit` must be of length N.")
        else:
            if len(m_unit) != N:
                raise ValueError("`m_unit` must be of length N.")

    def m_each(self, N, blocks, clusters, num_arms, conditions, simple, m_each):
        if any(m_each < 0):
            raise ValueError(
                "The number of units assigned to all conditions must be nonnegative."
            )
        if sum(m_each) != N:
            raise ValueError(
                "The sum of the number assigned to each condition (m_each) must equal the total number of units (N)."
            )

    def block_prob(self, N, blocks, clusters, num_arms, conditions, simple, block_prob):
        if any(block_prob_each < 0) or any(block_prob_each > 1):
            raise ValueError(
                "The probabilities of assignment to treatment must be between 0 and 1 for all blocks."
            )
        if len(block_prob_each) != blocks.shape[0]:
            raise ValueError(
                "If specified, block_prob_each should have the same number of rows as there are unique blocks in blocks."
            )
        if clusters is None and any(
            np.apply_along_axis(
                lambda x: np.isclose(np.sum(x), 1) != True, 1, block_prob_each
            )
        ):
            raise ValueError("If specified, each row of block_prob_each must sum to 1.")

    def block_m(self, block_m):
        if len(block_m) != blocks.attr("N_blocks"):
            raise ValueError(
                "If specified, block_m should have the same length as there are unique blocks in blocks."
            )

        if any(block_m > blocks.attr("N_per_block")) or any(block_m < 0):
            raise ValueError(
                "The number of units assigned to treatment within a block must be nonnegative and not exceed the total number units within the block."
            )

    def block_prob_each(
        self, N, blocks, clusters, num_arms, conditions, simple, block_prob_each
    ):
        if block_m_each.shape[0] != blocks.attr("N_blocks"):
            raise ValueError(
                "If specified, block_m_each should have the same number of rows as there are unique blocks in blocks."
            )
        if clusters is None and any(
            np.sum(block_m_each, axis=1) != blocks.attr("N_per_block")
        ):
            raise ValueError(
                "If specified, each row of block_m_each must sum to the number of units in the corresponding block."
            )
