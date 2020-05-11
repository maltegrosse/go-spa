Go-Solar Position Algorithm (SPA)
=======================================
[![Go Report Card](https://goreportcard.com/badge/github.com/maltegrosse/go-spa)](https://goreportcard.com/report/github.com/maltegrosse/go-spa)
[![GoDoc](https://godoc.org/github.com/maltegrosse/go-spa?status.svg)](https://pkg.go.dev/github.com/maltegrosse/go-spa)
![Go](https://github.com/maltegrosse/go-spa/workflows/Go/badge.svg) 

NREL's Solar Position Algorithm (SPA) calculates the solar zenith and azimuth angles in the period from the year -2000 to 6000, with uncertainties of +/- 0.0003 degrees based on the date, time, and location on Earth. 

(Reference: Reda, I.; Andreas, A., Solar Position Algorithm for Solar Radiation Applications, Solar Energy. Vol. 76(5), 2004; pp. 577-589). 
## Installation

This packages requires Go 1.13. If you installed it and set up your GOPATH, just run:

`go get -u github.com/maltegrosse/go-spa`

## Usage

You can find some examples in the [examples](examples) directory.

Please visit https://midcdmz.nrel.gov/spa/ for additional information.

Some additional helper functions have been added to the original application logic.
## Notes


|       | NREL spa_tester.c   | SPA GO    |  
|---------------|-------|-------|
| Julian Day     | 2452930.312847  | 2452930.312847  | 
| L          | 24.01826  | 24.018262  | 
| B      | -0.0001011219  | -0.000101121925  | 
| R            |  0.996542  | 0.996542  | 
| H  |  11.105902  |  11.105902  | 
| Delta Psi          | -0.003998404  | -0.003998404303  | 
| Delta Epsilon         | 0.001666568  | 0.001666568177  | 
| Epsilon          | 23.440465 | 23.440465  | 
| Zenith     | 50.111622 | 50.111622 | 
| Azimuth     | 194.340241  | 194.340241 | 
| Incidence      | 25.187000  | 25.187000  | 
| Sunrise         | 2003-10-17 06:12:43  | 2003-10-17 06:12:43 | 
| Sunset           | 2003-10-17 17:20:19  | 2003-10-17 17:20:19 |



## License
**[NREL SPA License](https://midcdmz.nrel.gov/spa/#license)**

Adoption in Golang under **[MIT license](http://opensource.org/licenses/mit-license.php)** 2020 © Malte Grosse.

