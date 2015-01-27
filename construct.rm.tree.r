library(RSQLite)

pdb <- dbConnect(SQLite(), "ancre.sqlite")

patientSelect <- dbGetQuery(pdb,
	"SELECT p.PatientID, count(s.Accession) FROM Patients AS p

INNER JOIN Sequences as s
ON s.PatientID = p.PatientID

INNER JOIN Products as pr
ON pr.Accession = s.Accession

WHERE p.StartDate IS NOT NULL and p.SIV = 0 and p.Treated = 0 and p.SuperInfection = 0 and pr.Product = 'env'
GROUP BY p.PatientID")

# Print a list of patients with the env gene sampled
cat(Reduce(paste0, apply(patientSelect, 1, function(s) sprintf("Patient %d has %d env sequences\n", s[1],s[2]))))

patient <- -1
while(! (as.numeric(patient) %in% patientSelect$PatientID))
  patient <- readline("enter a patient's number: ")

# 

tmpQuery <- sprintf("SELECT pr.Seq, s.Days FROM Patients AS p
INNER JOIN Sequences as s
ON s.PatientID = p.PatientID

INNER JOIN Products as pr
ON pr.Accession = s.Accession

WHERE p.StartDate IS NOT NULL and p.SIV = 0 and p.Treated = 0 and p.SuperInfection = 0 and pr.Product = 'env' and p.PatientID = '%d' ", as.numeric(patient))

sequences <- dbGetQuery(pdb, tmpQuery)
