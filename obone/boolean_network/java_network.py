import subprocess


def stepminer(
    bv_file: str,
    file_rl: str = "filler.rl",
    p_thr: float = 0.1,
    s_thr: float = 3.0,
    d_thr: float = 0.05,
) -> str:
    subprocess.run(
        [
            "java",
            "-cp",
            "stepminer-1.1.jar",
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrix",
            file_rl,
            bv_file,
            "false_filler.ph",
            "All",
            p_thr,
            s_thr,
            d_thr,
        ]
    )
    subprocess.run(
        [
            "java",
            "-cp",
            "stepminer-1.1.jar",
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrixFill",
            file_rl,
        ]
    )
    subprocess.run(
        [
            "java",
            "-cp",
            "stepminer-1.1.jar",
            "-Xms64m",
            "-Xmx10G",
            "tools.CustomAnalysis",
            "boolean",
            "bitMatrixFillStats",
            file_rl,
        ]
    )
    return file_rl
