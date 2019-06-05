#Data cleaning and prep for Sophie

# load data
BefSchnitt <- read.csv2('source_data/Schnittbefunde4.csv', sep = "\t", dec = ",")

# bringing the periods in the right order
BefSchnitt$epoche <- ordered(BefSchnitt$epoche, levels = c("Neolithikum", "Bronzezeit", "Eisenzeit", "undat."))
BefSchnitt$abschnitt <- ordered(BefSchnitt$abschnitt, levels =c('aelter', 'frueh','mittel', 'spaet'))
BefSchnitt$kultur <- ordered(BefSchnitt$kultur, levels=c("Baalberge", "Bernburg", "Kugelamphorenkultur", "Trichterbecher", "Salzmuende", "Schnurkeramik", "Glockenbecher", "Aunjetitz", "Saalemuendungsgruppe", "sBZ/aeEZ", "Hausurnenkultur"))

# subset settlement features
siedlung_BefSchnitt <- subset(BefSchnitt, BefSchnitt$befundart=='Siedlung')

#Subset the settlement features for relevant cultures
# there are 389 features dated as late bronze age and early iron age, that are being analysed in both epochs
siedlung_SnK_Bef <- subset(siedlung_BefSchnitt, siedlung_BefSchnitt$kultur=='Schnurkeramik')
siedlung_SnK_Bef <-siedlung_SnK_Bef[c('PlcIndex', 'Easting', 'Northing')]
write.csv(siedlung_SnK_Bef, "/home/sophie/Dokumente/Konferenzen/percolatransect/source_data/BefundesNL-minimal.csv", row.names = FALSE)

siedlung_fBZ_Bef <- subset(siedlung_BefSchnitt, siedlung_BefSchnitt$epoche=='Bronzezeit' & siedlung_BefSchnitt$abschnitt=='frueh')
siedlung_fBZ_Bef <-siedlung_fBZ_Bef[c("PlcIndex", 'Easting', 'Northing')]
write.csv(siedlung_fBZ_Bef, "/home/sophie/Dokumente/Konferenzen/percolatransect/source_data/BefundefBZ-minimal.csv")

siedlung_sBZ_Bef <-subset(siedlung_BefSchnitt, siedlung_BefSchnitt$epoche=='Bronzezeit' & siedlung_BefSchnitt$abschnitt=='spaet')
siedlung_sBZ_Bef <-siedlung_sBZ_Bef[c('PlcIndex', 'Easting', 'Northing')]
write.csv(siedlung_sBZ_Bef, "/home/sophie/Dokumente/Konferenzen/percolatransect/source_data/BefundesBZ-minimal.csv")


#siedlung_sBZ_aeEZ_Bef <- subset(siedlung_BefSchnitt, siedlung_BefSchnitt$epoche=='Eisenzeit' & siedlung_BefSchnitt$abschnitt=='ael
