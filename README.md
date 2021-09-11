# FEDDLib

Using AceGen element in the FEDDLib.
Root of modifications is the unsteadyTPM problem.
A copy is created in the ```feddlib/problems/Ace``` directoy for further modifications.
A copy of the specific problem "TPM" is created ```feddlib/problems/specific/ACE*``` for further modifications. The new specific problem class is called 
```
class ACE : public Problem<SC,LO,GO,NO>
```
Within the ```feddlib/core/FE``` structure, the ```void FE<SC,LO,GO,NO>::assemblyAceGenTPM``` was copied into the file ```feddlib/core/FE/Ace_def.hpp``` and renamed to ```void FE<SC,LO,GO,NO>::assemblyAceGenACE``` for further modifications. The file is simply *#included* into the ```FE_def.hpp```.

The AceLayer project is implemeted *dirty* for now.
The file ```feddlib/core/FE/ace_layer/ace_layer.hpp``` contains paths which are
system specific and also need to be constumized in install case.

## Notes:
 * AceGenElement class now available on every rank!