package org.scalanlp.radical

/**
  * Created by oar on 6/22/16.
  */
object Timer {

    var timeAtStart:Double = 0.0
    var timeAtStop:Double = 0.0

    def start: Unit = { timeAtStart = System.currentTimeMillis() }
    def stop:  Unit = { timeAtStop = System.currentTimeMillis() }
    def report: String = "Time: "+(timeAtStop-timeAtStart)+" millis."
    def printReport = System.out.println(report)
}
