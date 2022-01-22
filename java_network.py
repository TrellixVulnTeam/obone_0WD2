from cgi import test
import subprocess
import sys

$bitvector
$FILE
$P_THR
$S_THR
$D_THR


subprocess.run(
    java -cp stepminer-1.1.jar -Xms64m -Xmx10G tools.CustomAnalysis \
        boolean bitMatrix $FILE.rl \
        $bitvector \
        false_filler.ph All $P_THR $S_THR $D_THR
)

subprocess.run(
    java -cp stepminer-1.1.jar -Xms64m -Xmx10G tools.CustomAnalysis \
    boolean bitMatrixFill $FILE.rl
)

subprocess.run(
    java -cp stepminer-1.1.jar -Xms64m -Xmx10G tools.CustomAnalysis \
    boolean bitMatrixFillStats $FILE.rl
)

GSM PROBEID    000111002220000

/booleanfs2/sahoo/Data/Brain/RNASeq