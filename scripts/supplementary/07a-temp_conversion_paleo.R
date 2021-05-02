# Temperature calculation as per Westerhold et al. (2020)

# Last 66 million years ---------------------------------------------------

westerhold2020 <- read.csv("data/Table_34_Westerhold2020.csv")

# Equations from Westerhold 2020
# Deep Sea Temperature Tdo 
# 0.000 to 3.660 Tdo(°C) = 1 –4.4 * ((δ18O (‰) –3.25) / 3)
# 3.600 to 34.025Tdo(°C) = 5 –8 * ((δ18O (‰) –1.75) / 3)
# 34.025 to 67.000Tdo(°C) = (−4 * δ18O (‰)) + 12
# 
# Surface air temperature change TS 
# 0.000 to 1.810 Ts(°C) = 2 * Tdo + 12.25
# 1.810 to 5.330 Ts(°C) = 2.5 * Tdo + 12.15
# 5.330 to 67.000 Ts(°C) = Tdo + 14.15
# 
# Temperature anomaly with respect to average global temperature from 1961-1990
# Delta Temperature = Surface air temperature –14.15(Holocene mean temperature)

westerhold2020$d18O <- westerhold2020$ISOBENd18oLOESSsmoothLongTerm                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

#calculate deep sea temp

#+assign categories
westerhold2020$cat <- 1
westerhold2020$cat[westerhold2020$age_tuned > 3.660] <- 2
westerhold2020$cat[westerhold2020$age_tuned > 34.025] <- 3

westerhold2020$Tdo <- NA

westerhold2020$Tdo[westerhold2020$cat == 1] <- 1 -4.4 * ((westerhold2020$d18O[westerhold2020$cat == 1] -3.25) / 3)
westerhold2020$Tdo[westerhold2020$cat == 2] <- 5 -8 * ((westerhold2020$d18O[westerhold2020$cat == 2] -1.75) / 3)
westerhold2020$Tdo[westerhold2020$cat == 3] <- (-4 * westerhold2020$d18O[westerhold2020$cat == 3]) + 12

# calculate sea air temperature

#assign categories
westerhold2020$cat <- 1
westerhold2020$cat[westerhold2020$age_tuned > 1.810] <- 2
westerhold2020$cat[westerhold2020$age_tuned > 5.330] <- 3

westerhold2020$Ts <- NA

westerhold2020$Ts[westerhold2020$cat == 1] <- 2 * westerhold2020$Tdo[westerhold2020$cat == 1] + 12.25
westerhold2020$Ts[westerhold2020$cat == 2] <- 2.5 * westerhold2020$Tdo[westerhold2020$cat == 2] + 12.15
westerhold2020$Ts[westerhold2020$cat == 3] <- westerhold2020$Tdo[westerhold2020$cat == 3] + 14.15

plot(westerhold2020$age_tuned, westerhold2020$Ts, type="l", xlim=c(67,0))

# output
df <- westerhold2020[,c("age_tuned", "Ts")]
colnames(df) <- c("time", "temp")

write.table(df, file="ooet_pyrate/DES_input_data/WesterholdTemp.txt", sep="\t", row.names = F)
