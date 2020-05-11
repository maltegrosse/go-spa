package spa

// SPAFunctions defines the SPA functionalities
type SPAFunctions uint32

// enumeration for function codes to select desired final outputs from SPA
//go:generate stringer -type=SPAFunctions
const (
	SpaZa    SPAFunctions = 0 //calculate zenith and azimuth
	SpaZaInc SPAFunctions = 1 //calculate zenith, azimuth, and incidence
	SpaZaRts SPAFunctions = 2 //calculate zenith, azimuth, and sun rise/transit/set values
	SpaAll   SPAFunctions = 3 //calculate all SPA output values
)
