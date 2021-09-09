#include <ace_layer.hpp>
#include <iostream>

int AceFEMDataLayer::time_step(double CurrentTimeIncrement, double CurrentLoadMultiplier)
{
    // update step dependent integer data
    this->SystemIntegerData[ID_Step] += 1;
    this->SystemIntegerData[ID_Iteration] = 0;

    // update step dependent real data
    this->SystemRealData[RD_TimeIncrement] = CurrentTimeIncrement;
    this->SystemRealData[RD_Time] += CurrentTimeIncrement;
    this->SystemRealData[RD_Multiplier] = CurrentLoadMultiplier;

    return 0;
};