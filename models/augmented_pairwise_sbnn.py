"""Selected predictor-augmented SBNN used for full-background transfer."""

from augmented_sbnn_varpb_experiment import Config, fit_pairwise


MODEL_NAME = "augmented_pairwise_sbnn"
INPUT_FEATURES = 54
HIDDEN = (200, 100)
RANKS = (50, 50, 1)
KL_WEIGHT = 0.001
POSTERIOR_SCALE_MULTIPLIER = 0.1
PRIOR_SCALE = 1.0
POSTERIOR_DRAWS = 128


def fit_augmented_pairwise_sbnn(
    x_train,
    y_train,
    x_validation,
    y_validation,
    *,
    seed,
    epochs=30,
    pairs_per_epoch=40000,
):
    """Fit the configuration selected on the gene-disjoint validation split."""
    config = Config(
        RANKS,
        KL_WEIGHT,
        POSTERIOR_SCALE_MULTIPLIER,
        PRIOR_SCALE,
    )
    return fit_pairwise(
        x_train,
        y_train,
        x_validation,
        y_validation,
        config,
        seed=seed,
        epochs=epochs,
        pairs_per_epoch=pairs_per_epoch,
    )


def predict_augmented_pairwise_sbnn(model, x, *, draws=POSTERIOR_DRAWS, seed):
    """Return posterior mean, mutual information, and predictive deviation."""
    return model.predict_statistics(x, draws=draws, seed=seed)
