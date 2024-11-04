
from cc3d import CompuCellSetup
        

from ClusterCountingAlgorithmSteppables import ClusterCountingAlgorithmSteppable

CompuCellSetup.register_steppable(steppable=ClusterCountingAlgorithmSteppable(frequency=1))


CompuCellSetup.run()
