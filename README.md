# MAESPA-growth

This is a modified version of the model MAESPA that implements a simple growth module based on the growth-maintenance respiration approach, partitioning coefficients that change during the season,  a simple phenological approach based on thermal time and events describing harvest and pruning.

The implementation of these modules is achieved by modifying the daily loop of MAESPA and converting the maespa.f90 file into a subroutine callable from maespa_growth.f90 which now takes care of looping through the days. Note that the model only works for single-species, regular stands where all trees are assumed to have the same dimensions and properties.

# TO DO

- Pruning of root system when aboveground part is removed (to maintain root:shoot balance)
- Transfer 2-year-old shoots to wood (i.e. when leaves are shed the shoots become part of structural material)
