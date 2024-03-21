import coinstac
import run_vbm as vbm

# Start the computation, since this is preprocessing we can just pass the same script twice
# this should probably be cleaned up in the future so thats not necessary
coinstac.start(vbm.start, vbm.start)