import java.util.Scanner

case class EquationSystem(X: Vector[Vector[Double]], Y: Vector[Double]) {
  override def toString: String = X.indices.map(i => X(i).map(_.floor).mkString(" ") + "  " + Y(i).floor).mkString("\n")

  def updatedY(idx: Int, value: Double) = EquationSystem(X, Y.updated(idx, value))

  def size: Int = Y.length

  def showY: Unit = println(Y.map(n => "%.5f".format(n)).mkString(" "))
}

case class Row(coefficients: Vector[Double], Y: Double)


object GaussianElimination {

  /**
   * Input looks like:
   * 2
   * 1 1 3
   * 2 3 7
   */
  def main(args: Array[String]): Unit = {
    val equationSystem = scan

    def solve: EquationSystem => EquationSystem = eliminate _ andThen doBackPropagation

    val solved = solve(equationSystem)
    solved.showY
  }

  def eliminate(system: EquationSystem): EquationSystem = (0 until system.size - 1).foldLeft(system)((equation, i) => {
    val mr = maxRow(equation, i, i)
    val swapped = swap(equation, (i, mr))
    (i + 1 until system.size).foldLeft(swapped)((equation, j) => updateRow(equation, i, j))
  })

  def swap(equation: EquationSystem, swapPos: (Int, Int)): EquationSystem = {
    val equationUpdated = equation.X
      .updated(swapPos._1, equation.X(swapPos._2))
      .updated(swapPos._2, equation.X(swapPos._1))

    val valuesUpdated = equation.Y
      .updated(swapPos._1, equation.Y(swapPos._2))
      .updated(swapPos._2, equation.Y(swapPos._1))

    EquationSystem(equationUpdated, valuesUpdated)
  }

  def maxRow(equation: EquationSystem, startRow: Int, column: Int): Int = (startRow until equation.size).foldLeft(startRow)(
    (maxRow, row) => if (equation.X(row)(column).abs > equation.X(maxRow)(column).abs) row else maxRow
  )

  def updateRow(e: EquationSystem, sourceRow: Int, targetRow: Int): EquationSystem = {
    val rowToSource = e.X(sourceRow)
    val rowToUpdate = e.X(targetRow)
    val column = sourceRow //column number same as source number since we tackling diagonal of matrix

    val factor = rowToUpdate(column) / rowToSource(column)

    val updatedRow =
      rowToUpdate
        .updated(column, 0.0)
        .zip(rowToUpdate.indices.toVector)
        .map { case (value, idx) => if (idx > column && factor.abs > 0.00001) value - factor * rowToSource(idx) else value }

    log(s"${rowToUpdate.mkString(" ")} updated to ${updatedRow.mkString(" ")} with ${"%.3f".format(factor)} factor")

    val updatedCoefficients = e.X.updated(targetRow, updatedRow)
    val updatedValues = e.Y.updated(targetRow, e.Y(targetRow) - factor * e.Y(sourceRow))

    EquationSystem(updatedCoefficients, updatedValues)
  }

  def log(s: => String) = {
        println(s)
  }

  def scan: EquationSystem = {
    val scanner = new Scanner(System.in)
    val size = scanner.nextInt

    val X = Array.ofDim[Double](size, size)
    val Y = new Array[Double](size)

    for {
      row <- 0 until size
      col <- 0 to size
    } {
      if (col != size) X(row)(col) = scanner.nextInt
      else Y(row) = scanner.nextInt
    }

    EquationSystem(X.map(_.toVector).toVector, Y.toVector)
  }

  def doBackPropagation(system: EquationSystem): EquationSystem = system.Y.indices.reverse.foldLeft(system)((solution, idx) => {
    val sum = (idx + 1 until solution.size).map(j => solution.X(idx)(j) * solution.Y(j)).sum
    val varValue = (solution.Y(idx) - sum) / solution.X(idx)(idx)
    solution.updatedY(idx, varValue)
  })


}
