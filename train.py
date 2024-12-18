import polars as pl
import argparse


def train(fname: str):
    df = pl.read_csv(fname, separator="\t")

    df = (
        df.group_by(["ref_base_ctx_enc", "base_enc", "op"])
        .agg([pl.col("cnt").sum()])
        .with_columns(
            [
                (
                    pl.col("cnt")
                    / pl.col("cnt").sum().over([pl.col("ref_base_ctx_enc")])
                ).alias("freq")
            ]
        )
    )

    df = (
        df.with_columns(
            [
                pl.when(pl.col("op").eq("I"))
                .then(pl.format("+{}", pl.col("base_enc")))
                .otherwise(pl.col("base_enc"))
                .alias("base_enc")
            ]
        )
        .with_columns([pl.col("base_enc").fill_null("")])
        .select(
            [
                pl.format("{}->{}={}", "ref_base_ctx_enc", "base_enc", "freq").alias(
                    "##res"
                )
            ]
        )
    )

    df.write_csv("model.param")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("train_ins")
    train(parser.parse_args().train_ins)
