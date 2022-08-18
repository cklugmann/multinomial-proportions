import os
import matplotlib.pyplot as plt
from multinomial_estimation import MultinomialEstimator


def main():

    out_dir = "out"
    os.makedirs(out_dir, exist_ok=True)

    MultinomialEstimator.plot_sample_size(
        num_cells=3,
        population_size=10000
    )

    fig = plt.gcf()
    fig.savefig(os.path.join(out_dir, "sample_size.png"), bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
