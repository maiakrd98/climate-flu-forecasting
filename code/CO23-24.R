# reset weeks to start at 40
week = (ILINet$WEEK-40)%%52

# ILI times percent positive for 2023-2024
flu23 = ILINet$X.UNWEIGHTED.ILI*WHO_NREVSS_Clinical_Labs$PERCENT.POSITIVE

plot(week,log(flu23))

# ILI times percent positive for 2022-2023
flu22 = ILINet22$X.UNWEIGHTED.ILI*WHO_NREVSS_Clinical_Labs_22$PERCENT.POSITIVE

plot(week,log(flu22))

# ILI times percent positive for 2021-2022
flu21 = ILINet21$X.UNWEIGHTED.ILI*WHO_NREVSS_Clinical_Labs_21$PERCENT.POSITIVE

plot(week,log(flu21))

