package org.scalanlp.radical



/**
  * Created by oar on 7/2/16.
  * Ordered grid with binary search for the insertion index.
  */
final class DoubleGrid(val grid:Vector[Double]) {

    val lastIdx = grid.length-1
    def length = grid.length
    /** grid point i. */
    def apply(i:Int) = grid(i)

    def gridIsSorted =
        (0 until lastIdx).map(i => grid(i) < grid(i+1)).reduce(_&_)
    assert(gridIsSorted,"Grid values not sorted:\n"+grid+"\n")

    /* @return index j such that grid(j)<=u<grid(j+1).
     * if u<grid[0] returns -1
     * if u>=grid[lastIdx] returns lastIdx-1
     *
     * If the grid values are used as bin boundaries in a histogram, this is the number of the bin
     * the value u lands in (except if u<=grid(0)).
     */
    def leftIndex(u:Double):Int =

        if (u<grid(0)) -1 else if (u>=grid(lastIdx)) lastIdx-1
        else {

            var a=0; var b=lastIdx
            var m:Int = (a+b)/2
            // maintain  state grid(a)<=u<grid(b)
            while(b-a>1){
                if (grid(m)<=u) a=m else b=m
                m = (a+b)/2
            }
            a
        }

}
object DoubleGrid {

    /** factory method. */
    def apply(grid:Vector[Double]) = new DoubleGrid(grid)
}