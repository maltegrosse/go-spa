package main

import (
	"fmt"
	"github.com/maltegrosse/go-spa"
	"os"
	"text/tabwriter"
	"time"
)

func main() {

	deltaUt1 := 0.
	deltaT := 67.
	longitude := -105.1786
	latitude := 39.742476
	elevation := 1830.14
	pressure := 820.
	temperature := 11.
	slope := 30.
	azmRotation := -10.
	atmosRefract := 0.5667

	loc, err := time.LoadLocation("America/Los_Angeles")
	if err != nil {
		fmt.Println(err)
		return
	}
	dt := time.Date(2003, 10, 17, 12, 30, 30, 0, loc)
	spa, err := spa.NewSpa(dt, latitude, longitude, elevation, pressure, temperature, deltaT, deltaUt1, slope, azmRotation, atmosRefract)
	if err != nil {
		fmt.Println(err)
		return
	}

	writer := tabwriter.NewWriter(os.Stdout, 0, 8, 1, '\t', tabwriter.AlignRight)
	_, err = fmt.Fprintln(writer, "- \tNREL spa_tester.c\tGO SPA")
	if err != nil {
		fmt.Println(err)
		return
	}

	_, err = fmt.Fprintln(writer, "Julian Day", "\t", "2452930.312847", "\t", fmt.Sprintf("%.6f", spa.GetJd()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "L", "\t", "24.01826", "\t", fmt.Sprintf("%.6f", spa.GetL()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "B", "\t", "-0.0001011219", "\t", fmt.Sprintf("%.12f", spa.GetB()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "R", "\t", "0.996542", "\t", fmt.Sprintf("%.6f", spa.GetR()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "H", "\t", "11.105902", "\t", fmt.Sprintf("%.6f", spa.GetH()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}

	_, err = fmt.Fprintln(writer, "Delta Psi", "\t", "-0.003998404", "\t", fmt.Sprintf("%.12f", spa.GetDelPsi()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "Delta Epsilon", "\t", "0.001666568", "\t", fmt.Sprintf("%.12f", spa.GetDelEpsilon()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "Epsilon", "\t", "23.440465", "\t", fmt.Sprintf("%.6f", spa.GetEpsilon()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "Zenith", "\t", "50.111622", "\t", fmt.Sprintf("%.6f", spa.GetZenith()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}

	_, err = fmt.Fprintln(writer, "Azimuth", "\t", "194.340241", "\t", fmt.Sprintf("%.6f", spa.GetAzimuth()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}
	_, err = fmt.Fprintln(writer, "Incidence", "\t", "25.187000", "\t", fmt.Sprintf("%.6f", spa.GetIncidence()), "\t")
	if err != nil {
		fmt.Println(err)
		return
	}

	err = writer.Flush()
	if err != nil {
		fmt.Println(err)
		return
	}

	// use spa.Calculate() after the initial input values changed...

	fmt.Println(spa.GetSunrise())
	fmt.Println(spa.GetSunset())

	/////////////////////////////////////////////
	// The output of this program should be:
	//
	//Julian Day:    2452930.312847
	//L:             2.401826e+01 degrees
	//B:             -1.011219e-04 degrees
	//R:             0.996542 AU
	//H:             11.105902 degrees
	//Delta Psi:     -3.998404e-03 degrees
	//Delta Epsilon: 1.666568e-03 degrees
	//Epsilon:       23.440465 degrees
	//Zenith:        50.111622 degrees
	//Azimuth:       194.340241 degrees
	//Incidence:     25.187000 degrees
	//Sunrise:       06:12:43 Local Time
	//Sunset:        17:20:19 Local Time
	//
	/////////////////////////////////////////////
}
