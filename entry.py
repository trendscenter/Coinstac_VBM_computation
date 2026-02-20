import coinstac
import run_vbm as vbm

# Start preprocessing callbacks for both local and remote phases.
coinstac.start(vbm.start, vbm.start)
