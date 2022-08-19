using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GemBox.Spreadsheet;
using System.IO;
namespace HW6_309704035
{
    class Program
    {
        static void Main()
        {               
                #region//////User input
            #region // Matrix XB
            int RowXB, ColXB;
            Console.Write("Input the row of the XB matrix: ");
            RowXB = Convert.ToInt32(Console.ReadLine());
            //Console.Write("Input the column of the XB matrix: ");
            //ColXB = Convert.ToInt32(Console.ReadLine());
            ColXB = 1;
            double[][] XB = new double[RowXB][];
            for (int i = 0; i < RowXB; i++)
            {
                XB[i] = new double[ColXB];

            }
            for (int i = 0; i < RowXB; i++)
            {
                for (int j = 0; j < ColXB; j++)
                {
                    double input;
                    Console.Write("Enter value for ({0},{1}): ", i, j);
                    while (!double.TryParse(Console.ReadLine(), out input))
                    {
                        Console.Write("Enter correct value for ({0},{1}): ", i, j);
                    }
                    XB[i][j] = input;
                }
            }
            Console.Write("\nThe XB matrix is :\n");
            PrintMatrix(XB);
            Console.WriteLine();
            #endregion
            #region // Matrix B
            int RowB, ColB;
            //Console.Write("Input the row of the B matrix: ");
            //RowB = Convert.ToInt32(Console.ReadLine());
            RowB = XB.Length;
            ColB = XB.Length;
            Console.WriteLine("Input the value of the B matrix: ");
            //Console.Write("Input the column of the B matrix: ");
            //ColB = Convert.ToInt32(Console.ReadLine());
            double[][] B = new double[RowB][];
            for (int i = 0; i < RowB; i++)
            {
                B[i] = new double[ColB];
            }
            for (int i = 0; i < RowB; i++)
            {
                for (int j = 0; j < ColB; j++)
                {
                    double input;
                    Console.Write("Enter value for ({0},{1}): ", i, j);
                    while (!double.TryParse(Console.ReadLine(), out input))
                    {
                        Console.Write("Enter correct value for ({0},{1}): ", i, j);
                    }
                    B[i][j] = input;
                }
            }
            Console.Write("\nThe B matrix is :\n");
            PrintMatrix(B);
            #endregion
            #region // Matrix B^-1
            Console.Write("\nThe B^-1 matrix is :\n");
            double[][] InverB = Inverse(B);
            PrintMatrix(InverB);
            Console.WriteLine();

            #endregion
            #region // Matrix N
            int RowN, ColN;
            Console.WriteLine();
            //Console.Write("Input the row of the N matrix: ");
            //RowN = Convert.ToInt32(Console.ReadLine());
            RowN = XB.Length;
            Console.Write("Input the column of the N matrix: ");
            ColN = Convert.ToInt32(Console.ReadLine());
            double[][] N = new double[RowN][];
            for (int i = 0; i < RowN; i++)
            {
                N[i] = new double[ColN];
            }
            for (int i = 0; i < RowN; i++)
            {
                for (int j = 0; j < ColN; j++)
                {
                    double input;
                    Console.Write("Enter value for ({0},{1}): ", i, j);
                    while (!double.TryParse(Console.ReadLine(), out input))
                    {
                        Console.Write("Enter correct value for ({0},{1}): ", i, j);
                    }
                    N[i][j] = input;
                }
            }
            Console.Write("\nThe N matrix is :\n");
            PrintMatrix(N);
            #endregion
            #region // Matrix Basic
            Console.WriteLine();
            Console.Write("Input value of the Basic index matrix(index of basic variables): ");
            Console.WriteLine();
            int[] Basic = new int[ColB];
            for (int i = 0; i < Basic.Length; i++)
            {
                int input;
                Console.Write("Enter value for Basic({0}): ", i);
                while (!int.TryParse(Console.ReadLine(), out input))
                {
                    Console.Write("Enter correct value for Basic index matrix: ", i);
                }
                Basic[i] = input;
            }
            Console.WriteLine();
            Console.Write("\nThe Basic index matrix is :\n");
            Printarray(Basic);
            #endregion
            #region // Matrix NonBasic
            Console.Write("Input value of the NonBasic index matrix(index of NonBasic variables): ");
            Console.WriteLine();
            int[] NonBasic = new int[ColN];
            for (int i = 0; i < NonBasic.Length; i++)
            {
                int input;
                Console.Write("Enter value for NonBasic({0}): ", i);
                while (!int.TryParse(Console.ReadLine(), out input))
                {
                    Console.Write("Enter correct value for NonBasic index matrix: ", i);
                }
                NonBasic[i] = input;
            }
            Console.WriteLine();
            Console.Write("\nThe NonBasic index matrix is :\n");
            Printarray(NonBasic);
            #endregion
            #region// Matrix c
            int Rowc;
            Console.WriteLine();
            //Console.Write("Input the row of the c matrix: ");
            //Rowc = Convert.ToInt32(Console.ReadLine());
            Console.Write("Input values of the c=[cB,cN] matrix: ");
            Console.WriteLine();
            Rowc = Basic.Length + NonBasic.Length;
            double[][] c = new double[Rowc][];
            for (int i = 0; i < Rowc; i++)
            {
                c[i] = new double[1];
            }
            for (int i = 0; i < Rowc; i++)
            {
                for (int j = 0; j < c[i].Length; j++)
                {
                    double input;
                    Console.Write("Enter value for ({0},{1}): ", i, j);
                    while (!double.TryParse(Console.ReadLine(), out input))
                    {
                        Console.Write("Enter correct value for ({0},{1}): ", i, j);
                    }
                    c[i][j] = input;
                }
            }
            Console.WriteLine();
            #endregion
                #endregion
                #region//測試懶人包
                //double[][] XB = new double[][]
                //    {
                //    new double[] { 1 },
                //    new double[] { 3 },
                //    new double[] { 5 }
                //    };
                //double[][] B = new double[][]
                //    {
                //    new double[] {1,0,0  },
                //    new double[] {0,1,0  },
                //    new double[] {0,0,1  }
                //    };
                //double[][] InverB = new double[][]
                //    {
                //    new double[] { 1,0,0 },
                //    new double[] { 0,1,0 },
                //    new double[] { 0,0,1 }
                //    };
                //double[][] N = new double[][]
                //    {
                //    new double[] { 1,-1 },
                //    new double[] { 2,-1 },
                //    new double[] { 0, 1 }
                //    };
                //double[][] c = new double[][] 
                //    {
                //    new double[] {4},
                //    new double[] {3},
                //    new double[] {0},
                //    new double[] {0},
                //    new double[] {0}
                //    };
                //int[] Basic = new int[] { 3, 4, 5 };
                //int[] NonBasic = new int[] { 1, 2 };

                //Console.ReadKey();
                #endregion
                #region//測試懶人包2
            //double[][] XB = new double[][]
            //    {
            //    new double[] { 5 },
            //    new double[] { 3 },
            //    };
            //double[][] B = new double[][]
            //    {
            //    new double[] {1,0  },
            //    new double[] {0,1  },
            //    };
            //double[][] InverB = new double[][]
            //    {
            //    new double[] { 1,0 },
            //    new double[] { 0,1 },
            //    };
            //double[][] N = new double[][]
            //    {
            //    new double[] { 2,1,1,3 },
            //    new double[] { 1,3,1,2 },
            //    };
            //double[][] c = new double[][]  //[cN cB]
            //    {
            //    new double[] {0},
            //    new double[] {0},
            //    new double[] {6},
            //    new double[] {8},
            //    new double[] {5},
            //    new double[] {9}
            //    };
            //int[] Basic = new int[] { 5,6 };
            //int[] NonBasic = new int[] { 1,2,3,4};

            //Console.ReadKey();
            #endregion
                #region// variables
                int pv = NonBasic.Length;
                int[] PrimalVariables = new int[pv];
                for (int i = 0; i < pv; i++)
                {
                    PrimalVariables[i] = i + 1;
                    //Console.Write("x" + PrimalVariables[i] + "\n");
                }
                int[] PrimalSlackVariables = new int[B.Length];
                for (int i = 0; i < B.Length; i++)
                {
                    PrimalSlackVariables[i] = pv + 1 + i;
                    //Console.Write("x"+PrimalSlackVariables[i] + "\n");
                }
                int[] DualVariables = new int[PrimalSlackVariables.Length];
                for (int i = 0; i < PrimalSlackVariables.Length; i++)
                {
                    DualVariables[i] = 1 + i;
                    //Console.Write("y" + DualVariables[i] + "\n");
                }
                int[] DualSlackVariables = new int[PrimalVariables.Length];
                for (int i = 0; i < PrimalVariables.Length; i++)
                {
                    DualSlackVariables[i] = PrimalSlackVariables.Length + 1 + i;
                    //Console.Write("y" + DualSlackVariables[i] + "\n");
                }
                #endregion
                #region// cB c N , cBT  cNT
            double[][] cB = new double[Basic.Length][];
                for (int i = 0; i < Basic.Length; i++)
                {
                    cB[i] = new double[1];
                }
                for (int i = 0; i < Basic.Length; i++)
                {
                    cB[i][0] = c[i][0];
                }
            double[][] cN = new double[NonBasic.Length][];
                for (int i = 0; i < NonBasic.Length; i++)
                {
                    cN[i] = new double[1];
                }
                for (int i = 0; i < NonBasic.Length; i++)
                {
                    cN[i][0] = c[i + Basic.Length][0];
                }
                PrintMatrix(cB);
                PrintMatrix(cN);
                #endregion      
                #region//A
                double[][]A =new double[B.Length][];
                for (int i=0; i < A.Length;i++ )
                {
                    A[i] = new double[B[i].Length+N[i].Length];
                }
                for (int i = 0; i <A.Length; i++)
                {
                    for (int j = 0; j < A[i].Length; j++)
                    {
                        if(j<B.Length)
                            A[i][j] = B[i][j];
                        if(j>=B.Length)
                            A[i][j] = N[i][j - B.Length];
                    }
                }
                //PrintMatrix(A);
                #endregion
                #region//matrix I
                double[][] I = new double[InverB.Length][];
                for (int i = 0; i < InverB.Length; i++)
                {
                    I[i] = new double[InverB[i].Length];
                }
                for (int i = 0; i < InverB.Length; i++)
                {
                    for (int j = 0; j < InverB[i].Length; j++)
                    {
                        if (i == j)
                            I[i][j] = 1;
                        else
                            I[i][j] = 0;
                    }
                }
                #endregion
                #region//ej entering ,ei leaving
                double[][] ej = new double[NonBasic.Length][];
                for (int r = 0; r < NonBasic.Length; r++)
                {
                    ej[r] = new double[1];
                }
                //ei leaving
                double[][] ei = new double[Basic.Length][];
                for (int r = 0; r < Basic.Length; r++)
                {
                    ei[r] = new double[1];
                }
                #endregion
                #region//Vk
                double[][] Vk = new double[1][];
                for (int r = 0; r < Vk.Length; r++)
                {
                    Vk[r] = new double[NonBasic.Length];
                }
                #endregion
                var InverBN = MatrixMult(InverB, N);
                //Console.ReadLine();
                //StreamWriter sw = new StreamWriter(@"D:\share\HW6_[309704035]\iteration.csv"); //C:\Users\kenny\Desktop\線性規劃\HW6_309704035\iteration.csv
                StreamWriter sw = new StreamWriter(@"C:\Users\kenny\Desktop\git-repos\Simple Pivot Tool for Linear Programming\iteration.csv");
                //StreamWriter的路徑請自己設
                int s;
                sw.WriteLine("steepest edge rule for Initial: min(Vk) ");
                for (s = 1; s < 99999999; s++)
                {
                    var zN = MatrixSubtraction(MatrixMult(MatrixTranspose(MatrixMult(InverB, N)), cB), cN);
                    var zNT = MatrixTranspose(zN);
                    #region//steepest_edge_rule
                    for (int col = 0; col < NonBasic.Length; col++)
                    {
                        //entering
                        Make_ej(ej, col);
                        //Delta_xB = InverB*N*ej
                        var steepest_Delta_xB = MatrixMult(MatrixMult(InverB, N), ej);
                        //CHOOSE LEAVING
                        double[] steepest_CL = new double[steepest_Delta_xB.Length];
                        for (int i = 0; i < steepest_Delta_xB.Length; i++)
                        {
                            steepest_CL[i] = steepest_Delta_xB[i][0] / XB[i][0];
                        }
                        Make_ei(ei, Array.IndexOf(steepest_CL, steepest_CL.Max()));
                        //aj=N*ej (ak)
                        var ak = MatrixColumnSelection(A, NonBasic[col]);
                        var akT = MatrixTranspose(ak);
                        //Delta_xi
                        var steepest_Delta_xi = steepest_Delta_xB[Array.IndexOf(steepest_CL, steepest_CL.Max())][0];                        
                        //InverE
                        var steepest_InverE = MatrixSubtraction(I, MatrixDivision(MatrixMult(MatrixSubtraction(steepest_Delta_xB, ei), MatrixTranspose(ei)), steepest_Delta_xi));
                        //PrintMatrix(steepest_InverE);
                        // new InverB
                        var steepest_InverB = MatrixMult(InverB, steepest_InverE);
                        Console.WriteLine();
                        // make Vk
                        var vv = MatrixMult(MatrixMult(akT, MatrixTranspose(steepest_InverB)), MatrixMult(steepest_InverB, ak));
                        if (zNT[0].All(r => r <= 0))
                        {
                            Vk[0][col] = vv[0][0];
                        }
                        if (zNT[0][col] > 0)
                        {
                            Vk[0][col] = 306783067830678; // 不用算它
                        }
                        if (s > 1)
                        {
                            sw.WriteLine("steepest edge rule for Iteration" + (s-1) + ": min(Vk) "); 
                        }
                        Console.WriteLine("V" + NonBasic[col] + "=" + Vk[0][col]);
                        sw.WriteLine("V" + NonBasic[col] + "=" + Vk[0][col] + ",");
                        Console.WriteLine();
                    }
                    #endregion
                    #region/////Entering variable aj matrix
                    int Enterj = NonBasic[Array.IndexOf(Vk[0], Vk[0].Min())];
                    int indexEnterj = Array.IndexOf(NonBasic, Enterj);
                    double[][] aj = new double[N.Length][];
                    for (int i = 0; i < N.Length; i++)
                    {
                        aj[i] = new double[] { N[i][indexEnterj] }; // row length
                    }
                    #endregion
                    #region//result of matrix delta_xB[InverB*aj=delta_xB]
                    double[][] delta_xBj = new double[InverB.Length][]; //C == m*n
                    for (int i = 0; i < InverB.Length; i++)
                    {
                        delta_xBj[i] = new double[aj.Rank];
                    }
                    delta_xBj = MatrixMult(InverB, aj);
                    #endregion
                    #region//CHOOSE LEAVING
                    double[] CL = new double[delta_xBj.Length];
                    for (int i = 0; i < delta_xBj.Length; i++)
                    {
                        CL[i] = delta_xBj[i][0] / XB[i][0];
                    }
                    double t = 1 / (CL.Max()); ///t = value of Entering variableXj e.g. X1
                    double k = CL.Max();
                    int indexlv = Array.IndexOf(CL, CL.Max());
                    int lv = Basic[indexlv];
                    #endregion
                    #region//Update basic variable solutions ~XB (current)
                    for (int i = 0; i < XB.Length; i++)
                    {
                        for (int j = 0; j < XB.Rank; j++)
                        {
                            XB[i][j] = XB[i][j] - t * delta_xBj[i][j];
                            //Console.Write(string.Format("{0} ", XB[i][j]) + "\n"); // 印出陣列的值
                        }
                    }
                    Console.WriteLine();
                    #endregion
                    #region//Update the index sets of basic and nonbasic variables:
                    Basic[indexlv] = Enterj; 
                    NonBasic[indexEnterj] = lv;
                    int[] DualBasic = new int[NonBasic.Length];
                    for (int i = 0; i < DualBasic.Length; i++)
                    {
                        if (NonBasic[i] < DualVariables.Length)
                        {
                            DualBasic[i] = NonBasic[i] + PrimalSlackVariables.Length;
                            //Console.Write( DualBasic[i] + "\n");
                        }
                        else
                        {
                            DualBasic[i] = NonBasic[i] - PrimalVariables.Length;
                            //Console.Write( DualBasic[i] + "\n");
                        }
                    }
                    int[] DualNonBasic = new int[Basic.Length];
                    for (int i = 0; i < DualNonBasic.Length; i++)
                    {
                        if (Basic[i] < DualVariables.Length)
                        {
                            DualNonBasic[i] = Basic[i] + PrimalSlackVariables.Length;
                            //Console.Write( DualNonBasic[i] + "\n");
                        }
                        else
                        {                                                          // x1 x2 x3 x4 x5 
                            DualNonBasic[i] = Basic[i] - PrimalVariables.Length;   // y4 y5 y1 y2 y3 
                            //Console.Write( DualNonBasic[i] + "\n");
                        }
                    }
                    //Console.Write("\nThe new Basic index matrix is :\n");
                    //Printarray(Basic);
                    //Console.Write("\nThe new NonBasic index matrix is :\n");
                    //Printarray(NonBasic);
                    #endregion
                    #region//Update the matrix of ~XB (ONE DIRECTION)
                    XB[indexlv][0] = t;       //Xj = t
                    #endregion
                    #region// cB c N , cBT  cNT
                    double temp = cB[indexlv][0];
                    cB[indexlv][0] = cN[indexEnterj][0];
                    cN[indexEnterj][0] = temp;
                    PrintMatrix(cB);
                    PrintMatrix(cN);
                    var cBT = MatrixTranspose(cB);
                    var cNT = MatrixTranspose(cN);
                    #endregion
                    #region//ei
                    //double[][] ei = new double[delta_xBj.Length][];
                    for (int i = 0; i < ei.Length; i++)
                    {
                        ei[i] = new double[delta_xBj[i].Length];
                        for (int j = 0; j < ei[i].Length; j++)
                        {
                            ei[i][j] = 0;
                        }
                    }
                    ei[indexlv][0] = 1;
                    //Console.Write("\nThe ei matrix is :\n");
                    //PrintMatrix(ei);
                    #endregion
                    #region//eiT
                    var eiT = MatrixTranspose(ei);
                    #endregion
                    #region// delta_xBi
                    double delta_xBi = delta_xBj[indexlv][0];
                    #endregion
                    #region caculate InverE matrix= I -(delta_xBj-ei)*eiT/delta_xi
                    //## caculate InverE matrix= I-(R1)*eiT/delta_xi=I-R2/delta_xi  
                    #endregion
                    #region// InverE  =I -(delta_xBj-ei)*eiT/delta_xBi
                    var InverE = MatrixSubtraction(I, MatrixDivision(MatrixMult(MatrixSubtraction(delta_xBj, ei), eiT), delta_xBi));
                    var E = MatrixTranspose(InverE);
                    //PrintMatrix(InverE);
                    #endregion
                    #region// newInverB = InverE*InverB
                    InverB = MatrixMult(InverE, InverB);
                    #endregion
                    #region// newN
                    for (int i = 0; i < N.Length; i++)
                    {
                        for (int j = 0; j < N[i].Length; j++)
                        {
                            if (j == indexEnterj)
                                N[i][indexEnterj] = B[i][indexlv];
                            else
                                N[i][j] = N[i][j];
                        }
                    }
                    #endregion
                    #region//cBT*newXb zN
                    zN = MatrixSubtraction(MatrixMult(MatrixTranspose(MatrixMult(InverB, N)), cB), cN);
                    zNT = MatrixTranspose(zN);
                    var Objvalue = MatrixMult(cBT, XB);
                    PrintMatrix(zN);
                    #endregion
                    #region// primal dictionary matrix
                    InverBN = MatrixMult(InverB, N);
                    Console.Write("\nThe primal dictionary matrix\t" + s + "\tis :\n\n");
                    Console.Write("\t");
                    double[][] DictionaryMatrix = new double[XB.Length + 1][];
                    for (int i = 0; i < XB.Length + 1; i++)
                    {
                        DictionaryMatrix[i] = new double[XB[0].Length + InverBN[0].Length];//XB[i].Length + InverBN[i].Length
                    }
                    DictionaryMatrix[0][0] = Objvalue[0][0];
                    for (int i = 1; i < DictionaryMatrix[0].Length; i++)
                    {
                        DictionaryMatrix[0][i] = -zNT[0][i - 1];
                    }
                    for (int i = 1; i < DictionaryMatrix.Length; i++)
                    {
                        DictionaryMatrix[i][0] = XB[i - 1][0];
                    }
                    for (int i = 1; i < DictionaryMatrix.Length; i++)
                    {
                        for (int j = 1; j < DictionaryMatrix[i].Length; j++)
                        {
                            DictionaryMatrix[i][j] = -InverBN[i - 1][j - 1];
                        }
                    }
                    for (int i = 0; i < NonBasic.Length; i++)
                    {
                        Console.Write("\tx" + NonBasic[i]);
                    }
                    Console.WriteLine();
                    for (int i = 0; i < DictionaryMatrix[0].Length; i++)
                    {
                        Console.Write("\t" + DictionaryMatrix[0][i]);
                    }
                    Console.WriteLine();
                    for (int i = 1; i < DictionaryMatrix.Length; i++)
                    {
                        Console.Write("x" + Basic[i - 1] + "\t");
                        for (int j = 0; j < DictionaryMatrix[i].Length; j++)
                        {
                            Console.Write(DictionaryMatrix[i][j] + "\t");
                        }
                        Console.WriteLine();
                    }
                    Console.Write("\nObjective value = " + Objvalue[0][0] + " ");
                    Console.Write(",Primal variables :");

                    for (int i = 0; i < PrimalVariables.Length; i++)
                    {
                        if (NonBasic.Any(r => r == PrimalVariables[i]))
                            Console.Write("x" + PrimalVariables[i] + " = " + 0 + "  ");  // B= 3,4,5 N= 1,2 XB = 1,3,5
                        if (Basic.Any(r => r == PrimalVariables[i]))
                            Console.Write("x" + PrimalVariables[i] + " = " + XB[Array.IndexOf(Basic, PrimalVariables[i])][0] + "  ");
                    }
                    Console.Write(",Primal slack variables :");
                    for (int i = 0; i < PrimalSlackVariables.Length; i++)
                    {
                        if (NonBasic.Any(r => r == PrimalSlackVariables[i]))
                            Console.Write("x" + PrimalSlackVariables[i] + " = " + 0 + "  ");
                        if (Basic.Any(r => r == PrimalSlackVariables[i]))
                            Console.Write("x" + PrimalSlackVariables[i] + " = " + XB[Array.IndexOf(Basic, PrimalSlackVariables[i])][0] + "  ");
                    }
                    Console.WriteLine();
                    #endregion
                    #region//dual dictionary matrix
                    var Dual = MatrixNegativeTranspose(DictionaryMatrix);
                    Console.Write("\nThe Dual Dictionary Matrix\t" + s + "\tis :\n\n");
                    //PrintMatrix(Dual);
                    Console.Write("\t");
                    for (int i = 0; i < DualNonBasic.Length; i++)
                    {
                        Console.Write("\ty" + DualNonBasic[i]);
                    }
                    Console.WriteLine();
                    for (int i = 0; i < Dual[0].Length; i++)
                    {
                        Console.Write("\t" + Dual[0][i]);
                    }
                    Console.WriteLine();
                    for (int i = 1; i < Dual.Length; i++)
                    {
                        Console.Write("y" + DualBasic[i - 1] + "\t");
                        for (int j = 0; j < Dual[i].Length; j++)
                        {
                            Console.Write(Dual[i][j] + "\t");
                        }
                        Console.WriteLine();
                    }
                    Console.Write("\nDual Objective value = " + -Objvalue[0][0] + " ");
                    Console.Write(",Dual variables :");
                    for (int i = 0; i < DualVariables.Length; i++)
                    {
                        if (DualNonBasic.Any(r => r == DualVariables[i]))
                            Console.Write("y" + DualVariables[i] + " = " + 0 + " ");
                        if (DualBasic.Any(r => r == DualVariables[i]))
                            Console.Write("y" + DualVariables[i] + " = " + Dual[Array.IndexOf(DualBasic, DualVariables[i]) + 1][0] + " ");
                    }
                    Console.Write(" ,Dual slack variables :");
                    for (int i = 0; i < DualSlackVariables.Length; i++)
                    {
                        if (DualNonBasic.Any(r => r == DualSlackVariables[i]))
                            Console.Write("y" + DualSlackVariables[i] + " = " + 0 + " ");
                        if (DualBasic.Any(r => r == DualSlackVariables[i]))
                            Console.Write("y" + DualSlackVariables[i] + " = " + Dual[Array.IndexOf(DualBasic, DualSlackVariables[i]) + 1][0] + " ");
                    }
                    Console.WriteLine();
                    sw.WriteLine();
                    #endregion
                    #region// csv file
                    //primal
                    sw.WriteLine("Primal Dictionary Matrix" + s + ":");
                    sw.Write("," + ",");
                    for (int i = 0; i < NonBasic.Length; i++)
                    {
                        sw.Write("x" + NonBasic[i] + ",");
                    }
                    sw.WriteLine();
                    for (int i = 0; i < DictionaryMatrix.Length; i++)
                    {
                        if (i == 0)
                            sw.Write(",");
                        else
                            sw.Write("x" + Basic[i - 1] + ",");
                        for (int j = 0; j < DictionaryMatrix[i].Length; j++)
                        {
                            sw.Write(DictionaryMatrix[i][j] + ",");
                        }
                        sw.WriteLine();
                    }
                    sw.Write("Objective value = " + Objvalue[0][0]);
                    sw.WriteLine();
                    sw.Write("prima variables:" + "," + ",");
                    for (int i = 0; i < PrimalVariables.Length; i++)
                    {
                        if (NonBasic.Any(r => r == PrimalVariables[i]))
                            sw.Write("x" + PrimalVariables[i] + " = " + 0 + ",");  // B= 3,4,5 N= 1,2 XB = 1,3,5
                        if (Basic.Any(r => r == PrimalVariables[i]))
                            sw.Write("x" + PrimalVariables[i] + " = " + XB[Array.IndexOf(Basic, PrimalVariables[i])][0] + ",");
                    }
                    sw.WriteLine();
                    sw.Write("primal slack variables :" + "," + ",");
                    for (int i = 0; i < PrimalSlackVariables.Length; i++)
                    {
                        if (NonBasic.Any(r => r == PrimalSlackVariables[i]))
                            sw.Write("x" + PrimalSlackVariables[i] + " = " + 0 + ",");
                        if (Basic.Any(r => r == PrimalSlackVariables[i]))
                            sw.Write("x" + PrimalSlackVariables[i] + " = " + XB[Array.IndexOf(Basic, PrimalSlackVariables[i])][0] + ",");
                    }
                    //dual
                    sw.WriteLine();
                    sw.WriteLine();
                    sw.WriteLine("Dual Dictionary Matrix" + s + ":");
                    sw.Write("," + ",");
                    for (int i = 0; i < DualNonBasic.Length; i++)
                    {
                        sw.Write("y" + DualNonBasic[i] + ",");
                    }
                    sw.WriteLine();
                    for (int i = 0; i < Dual.Length; i++)
                    {
                        if (i == 0)
                            sw.Write(",");
                        else
                            sw.Write("y" + DualBasic[i - 1] + ",");
                        for (int j = 0; j < Dual[i].Length; j++)
                        {
                            sw.Write(Dual[i][j] + ",");
                        }
                        sw.WriteLine();
                    }
                    sw.Write("Objective value = " + Objvalue[0][0]);
                    sw.WriteLine();
                    sw.Write("dual variables:" + "," + ",");
                    for (int i = 0; i < DualVariables.Length; i++)
                    {
                        if (DualNonBasic.Any(r => r == DualVariables[i]))
                            sw.Write("y" + DualVariables[i] + " = " + 0 + ",");
                        if (DualBasic.Any(r => r == DualVariables[i]))
                            sw.Write("y" + DualVariables[i] + " = " + Dual[Array.IndexOf(DualBasic, DualVariables[i]) + 1][0] + ",");
                    }
                    sw.WriteLine();
                    sw.Write("dual slack variables :" + "," + ",");
                    for (int i = 0; i < DualSlackVariables.Length; i++)
                    {
                        if (DualNonBasic.Any(r => r == DualSlackVariables[i]))
                            sw.Write("y" + DualSlackVariables[i] + " = " + 0 + ",");
                        if (DualBasic.Any(r => r == DualSlackVariables[i]))
                            sw.Write("y" + DualSlackVariables[i] + " = " + Dual[Array.IndexOf(DualBasic, DualSlackVariables[i]) + 1][0] + ",");
                    }
                    sw.WriteLine();
                    sw.WriteLine();
                    #endregion
                    if (zNT[0].All(r => r >= 0)) break;
                }
                sw.WriteLine();
                sw.WriteLine("Iteration " + s + " is Optimal");
                sw.Close();
                Console.WriteLine("\n\tOptimal");
                Console.ReadLine();   //要看console請打開
            }

        public static double[][] Inverse(double[][] Array)
        {
            int m = 0;
            int n = 0;
            m = Array.Length;
            n = Array[0].Length;
            double[][] array = new double[2 * m + 1][];
            for (int i = 0; i < 2 * m + 1; i++)
            {
                array[i] = new double[2 * n + 1];
            }
            for (int k = 0; k < 2 * m + 1; k++)  //初始化數組
            {
                for (int t = 0; t < 2 * n + 1; t++)
                {
                    array[k][t] = 0.00000000;
                }
            }
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    array[i][j] = Array[i][j];
                }
            }

            for (int k = 0; k < m; k++)
            {
                for (int t = n; t <= 2 * n; t++)
                {
                    if ((t - k) == m)
                    {
                        array[k][t] = 1.0;
                    }
                    else
                    {
                        array[k][t] = 0;
                    }
                }
            }
            //得到逆矩陣
            for (int k = 0; k < m; k++)
            {
                if (array[k][k] != 1)
                {
                    double bs = array[k][k];
                    array[k][k] = 1;
                    for (int p = k + 1; p < 2 * n; p++)
                    {
                        array[k][p] /= bs;
                    }
                }
                for (int q = 0; q < m; q++)
                {
                    if (q != k)
                    {
                        double bs = array[q][k];
                        for (int p = 0; p < 2 * n; p++)
                        {
                            array[q][p] -= bs * array[k][p];
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
            }
            double[][] NI = new double[m][];
            for (int i = 0; i < m; i++)
            {
                NI[i] = new double[n];
            }
            for (int x = 0; x < m; x++)
            {
                for (int y = n; y < 2 * n; y++)
                {
                    NI[x][y - n] = array[x][y];
                }
            }
            return NI;
        }//反矩陣
        private static void PrintMatrix(double[][] matrix)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    Console.Write(matrix[i][j] + "\t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }
        private static void WriteMatrix(double[][] matrix)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j < matrix[i].Length; j++)
                {
                    Console.Write(matrix[i][j] + "\n");
                }
            }
        }
        private static void Printarray(int[] array)
        {
            for (int i = 0; i < array.Length; i++)
            {
                Console.Write(array[i] + "\t");
            }
            Console.WriteLine();

        }
        private static void Printdoublearray(double[] doublearray)
        {
            for (int i = 0; i < doublearray.Length; i++)
            {
                Console.Write(doublearray[i] + "\t");
            }
            Console.WriteLine();

        }
        private static double[][] MatrixMult(double[][] matrix1, double[][] matrix2)
        {
            //matrix1是m*n矩陣，matrix2是n*p矩陣，則result是m*p矩陣
            int m = matrix1.Length, n = matrix2.Length, p = matrix2[0].Length;
            double[][] result = new double[m][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[p];
            }
            //矩陣乘法：c[i,j]=Sigma(k=1→n,a[i,k]*b[k,j])
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < p; j++)
                {
                    //對乘加法則
                    for (int k = 0; k < n; k++)
                    {
                        result[i][j] += (matrix1[i][k] * matrix2[k][j]);
                    }
                }
            }
            return result;
        }
        private static double[][] MatrixTranspose(double[][] matrix)
        {
            //矩陣中沒有元素的情況
            if (matrix.Length == 0)
            {
                return new double[][] { };
            }
            double[][] result = new double[matrix[0].Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix.Length];
            }
            //新矩陣生成規則： b[i,j]=a[j,i]
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result[0].Length; j++)
                {
                    result[i][j] = matrix[j][i];
                }
            }
            return result;
        }
        private static double[][] MatrixNegativeTranspose(double[][] matrix)
        {
            //矩陣中沒有元素的情況
            if (matrix.Length == 0)
            {
                return new double[][] { };
            }
            double[][] result = new double[matrix[0].Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix.Length];
            }
            //新矩陣生成規則： b[i,j]=a[j,i]
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result[0].Length; j++)
                {
                    result[i][j] = -matrix[j][i];
                }
            }
            return result;
        }
        private static double[][] MatrixSubtraction(double[][] matrix1, double[][] matrix2)
        {
            //生成一個與matrix1同型的空矩陣
            double[][] result = new double[matrix1.Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix1[i].Length];
            }
            //矩陣加法：把矩陣2各元素值加到矩陣1上，返回矩陣1
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result[i].Length; j++)
                {
                    result[i][j] = matrix1[i][j] - matrix2[i][j];
                }
            }
            return result;
        }
        private static double[][] MatrixDivision(double[][] matrix1, double delta)
        {
            double[][] result = new double[matrix1.Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix1[i].Length];
            }
            for (int i = 0; i < result.Length; i++)
            {
                for (int j = 0; j < result[i].Length; j++)
                {
                    result[i][j] = matrix1[i][j] / delta;
                }
            }
            return result;
        }
        private static double[][] MatrixColumnSelection(double[][] matrix1, int col)
        {
            double[][] result = new double[matrix1.Length][];
            for (int i = 0; i < matrix1.Length; i++)
            {
                result[i] = new double[1];
            }
            for (int i = 0; i < matrix1.Length; i++)
            {
                result[i][0] = matrix1[i][col-1];
            }
            return result;
        }
        private static double[][] MatrixRowSelection(double[][] matrix1, int Row)
        {
            double[][] result = new double[1][];
            for (int i = 0; i < matrix1.Length; i++)
            {
                result[i] = new double[matrix1[i].Length];
            }
            for (int j = 0; j < matrix1[j].Length; j++)
            {
                result[Row][j] = matrix1[Row][j];
            }
            return result;
        }
        private static double[][] Delta_xB(double[][] matrix, int col)
        {
            double[][] result = new double[matrix.Length][];
            for (int i = 0; i < matrix.Length; i++)
            {
                result[i] = new double[1];
            }
            for (int i = 0; i < matrix.Length; i++)
            {
                result[i][0] = matrix[i][col];
            }
            return result;
        }
        private static void Make_ej(double[][] matrix, int entering_index)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                matrix[i][0] = 0;
                matrix[entering_index][0] = 1;
            }
        }
        private static void Make_ei(double[][] matrix, int leaving_index)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                matrix[i][0] = 0;
                matrix[leaving_index][0] = 1;
            }
        }
    }
    }

