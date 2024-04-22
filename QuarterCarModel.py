#region imports
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import math
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg

#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
#endregion

#region class definitions
#region specialized graphic items
class MassBlock(qtw.QGraphicsItem):
    def __init__(self, CenterX, CenterY, width=30, height=10, parent=None, pen=None, brush=None, name='CarBody', mass=10):
        super().__init__(parent)
        self.x = CenterX
        self.y = CenterY
        self.pen = pen
        self.brush = brush
        self.width = width
        self.height = height
        self.top = self.y - self.height/2
        self.left = self.x - self.width/2
        self.rect = qtc.QRectF(self.left, self.top, self.width, self.height)
        self.name = name
        self.mass = mass
        self.transformation = qtg.QTransform()
        stTT = self.name +"\nx={:0.3f}, y={:0.3f}\nmass = {:0.3f}".format(self.x, self.y, self.mass)
        self.setToolTip(stTT)

    def boundingRect(self):
        bounding_rect = self.transformation.mapRect(self.rect)
        return bounding_rect

    def paint(self, painter, option, widget=None):
        self.transformation.reset()
        if self.pen is not None:
            painter.setPen(self.pen)  # Red color pen
        if self.brush is not None:
            painter.setBrush(self.brush)
        self.top = -self.height/2
        self.left = -self.width/2
        self.rect=qtc.QRectF( self.left, self.top, self.width, self.height)
        painter.drawRect(self.rect)
        self.transformation.translate(self.x, self.y)
        self.setTransform(self.transformation)
        self.transformation.reset()
        brPen=qtg.QPen()
        brPen.setWidth(0)
        painter.setPen(brPen)
        painter.setBrush(qtc.Qt.NoBrush)
        painter.drawRect(self.boundingRect())

class Wheel(qtw.QGraphicsItem):
    def __init__(self, CenterX, CenterY, radius=10, parent=None, pen=None, wheelBrush=None, massBrush=None, name='Wheel', mass=10):
        super().__init__(parent)
        self.x = CenterX
        self.y = CenterY
        self.pen = pen
        self.brush = wheelBrush
        self.radius = radius
        self.rect = qtc.QRectF(self.x - self.radius, self.y - self.radius, self.radius*2, self.radius*2)
        self.name = name
        self.mass = mass
        self.transformation = qtg.QTransform()
        stTT = self.name +"\nx={:0.3f}, y={:0.3f}\nmass = {:0.3f}".format(self.x, self.y, self.mass)
        self.setToolTip(stTT)
        self.massBlock = MassBlock(CenterX, CenterY, width=2*radius*0.85, height=radius/3, pen=pen, brush=massBrush, name="Wheel Mass", mass=mass)

    def boundingRect(self):
        bounding_rect = self.transformation.mapRect(self.rect)
        return bounding_rect
    def addToScene(self, scene):
        scene.addItem(self)
        scene.addItem(self.massBlock)

    def paint(self, painter, option, widget=None):
        self.transformation.reset()
        if self.pen is not None:
            painter.setPen(self.pen)  # Red color pen
        if self.brush is not None:
            painter.setBrush(self.brush)
        self.rect=qtc.QRectF(-self.radius, -self.radius, self.radius*2, self.radius*2)
        painter.drawEllipse(self.rect)
        self.transformation.translate(self.x, self.y)
        self.setTransform(self.transformation)
        self.transformation.reset()
        brPen=qtg.QPen()
        brPen.setWidth(0)
        painter.setPen(brPen)
        painter.setBrush(qtc.Qt.NoBrush)
        painter.drawRect(self.boundingRect())

#endregion

#region MVC for quarter car model
class CarModel():
    """
    I re-wrote the quarter car model as an object oriented program
    and used the MVC pattern.  This is the quarter car model.  It just
    stores information about the car and results of the ode calculation.
    """
    def __init__(self):
        """
        self.results to hold results of odeint solution
        self.t time vector for odeint and for plotting
        self.tramp is time required to climb the ramp
        self.angrad is the ramp angle in radians
        self.ymag is the ramp height in m
        """
        self.results = []
        self.tmax = 3.0  # limit of timespan for simulation in seconds
        self.t = np.linspace(0, self.tmax, 200)
        self.tramp = 1.0  # time to traverse the ramp in seconds
        self.angrad = 0.1
        self.ymag = 6.0 / (12 * 3.3)  # ramp height in meters.  default is 0.1515 m
        self.yangdeg = 45.0  # ramp angle in degrees.  default is 45
        self.results = None

        #set default values for the properties of the quarter car model
        self.m1 = 450 # mass of car body in kg
        self.m2 = 20  # mass of wheel in kg
        self.c1 = 4500  # damping coefficient in N*s/m
        self.k1 = 43421.46  # average spring constant for suspension in N/m
        self.k2 = 173685.825  # average spring constant for tire in N/m
        self.v = 120  # velocity of car in kph


        self.mink1 = 80000  #If I jack up my car and release the load on the spring, it extends about 3 inches
        self.maxk1 = 120000  #What would be a good value for a soft spring vs. a stiff spring?
        self.mink2 = 100000  #Same question for the shock absorber.
        self.maxk2 = 200000
        self.accel =None
        self.accelMax = 0
        self.accelLim = 2.0 * 9.81
        self.SSE = 0.0

class CarView():
    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        # unpack widgets with same names as they have on the GUI
        self.le_m1, self.le_v, self.le_k1, self.le_c1, self.le_m2, self.le_k2, self.le_ang, \
         self.le_tmax, self.chk_IncludeAccel = self.input_widgets

        self.gv_Schematic, self.chk_LogX, self.chk_LogY, self.chk_LogAccel, \
        self.chk_ShowAccel, self.lbl_MaxMinInfo, self.layout_horizontal_main = self.display_widgets

        # creating a canvas to draw a figure for the car model
        self.figure = Figure(tight_layout=True, frameon=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.layout_horizontal_main.addWidget(self.canvas)

        # axes for the plotting using view
        self.ax = self.figure.add_subplot()
        if self.ax is not None:
            self.ax1 = self.ax.twinx()

        self.buildScene()

    def updateView(self, model=None):
        self.le_m1.setText("{:0.2f}".format(model.m1))
        self.le_k1.setText("{:0.2f}".format(model.k1))
        self.le_c1.setText("{:0.2f}".format(model.c1))
        self.le_m2.setText("{:0.2f}".format(model.m2))
        self.le_k2.setText("{:0.2f}".format(model.k2))
        self.le_ang.setText("{:0.2f}".format(model.yangdeg))
        self.le_tmax.setText("{:0.2f}".format(model.tmax))
        stTmp="k1_min = {:0.2f}, k1_max = {:0.2f}\nk2_min = {:0.2f}, k2_max = {:0.2f}\n".format(model.mink1, model.maxk1, model.mink2, model.maxk2)
        stTmp+="SSE = {:0.2f}".format(model.SSE)
        self.lbl_MaxMinInfo.setText(stTmp)
        self.doPlot(model)

    def buildScene(self):
        # Create a scene object
        self.scene = qtw.QGraphicsScene()
        self.scene.setObjectName("MyScene")
        self.scene.setSceneRect(-200, -200, 400, 400)  # Set scene dimensions

        # Set the scene for the graphics view object
        self.gv_Schematic.setScene(self.scene)

        # Setup pens and brushes for drawing
        self.setupPensAndBrushes()

        # create and add the Car Body (Rectangle)
        car_body_rect = qtw.QGraphicsRectItem(-50, -100, 100, 20)
        car_body_rect.setPen(self.penWheel)
        car_body_rect.setBrush(self.brushMass)
        self.scene.addItem(car_body_rect)
        # mass label
        car_body_label = qtw.QGraphicsTextItem("mass = 450 kg", car_body_rect)
        car_body_label.setPos(-30, -120)
        self.scene.addItem(car_body_label)

        # add the wheel (Circle)
        wheel_circle = qtw.QGraphicsEllipseItem(-40, -50, 75, 75)
        wheel_circle.setPen(self.penWheel)
        wheel_circle.setBrush(self.brushWheel)
        self.scene.addItem(wheel_circle)
        # Add mass label to Wheel
        wheel_label = qtw.QGraphicsTextItem("mass = 20 kg", wheel_circle)
        wheel_label.setPos(-35, 20)  # Adjust the position based on your scene coordinates
        self.scene.addItem(wheel_label)

        # Create and add the Spring (Zigzag line)
        spring_lines = []
        spring_start_x, spring_start_y = -25, -100  # Starting point at the bottom center of the car body
        spring_end_x, spring_end_y = -25, -25  # Ending point at the top center of the wheel
        for i in range(5):  # Create 5 segments for the spring
            if i % 2 == 0:
                end_x = spring_start_x + 15  # Zig to the right
            else:
                end_x = spring_start_x - 15  # Zag to the left
            end_y = spring_start_y + (spring_end_y - spring_start_y) / 5  # Increment y position
            line = qtw.QGraphicsLineItem(spring_start_x, spring_start_y, end_x, end_y)
            line.setPen(self.penSpring)
            self.scene.addItem(line)
            spring_lines.append(line)
            spring_start_x, spring_start_y = end_x, end_y  # Update to new starting point for next segment

        # Create and add the Damper (Line with rectangle)
        damper_line = qtw.QGraphicsLineItem(25, -100, 25, -15)
        damper_line.setPen(self.penDamper)
        self.scene.addItem(damper_line)
        # Add rectangle to represent the damper housing
        damper_housing = qtw.QGraphicsRectItem(20, -75, 10, 5)
        damper_housing.setBrush(self.brushWheel)
        self.scene.addItem(damper_housing)


    def setupPensAndBrushes(self):
        # Define pens and brushes
        self.penWheel = qtg.QPen(qtg.QColor("orange"), 2)
        self.penSpring = qtg.QPen(qtg.QColor("black"), 2)
        self.penDamper = qtg.QPen(qtg.QColor("grey"), 2)
        self.brushWheel = qtg.QBrush(qtg.QColor.fromHsv(35,255,255, 64))
        self.brushMass = qtg.QBrush(qtg.QColor(200,200,200, 128))



    #def setupPensAndBrushes(self):
       # self.penWheel = qtg.QPen(qtg.QColor("orange"))
        #self.penWheel.setWidth(1)
       # self.brushWheel = qtg.QBrush(qtg.QColor.fromHsv(35,255,255, 64))
       # self.brushMass = qtg.QBrush(qtg.QColor(200,200,200, 128))

    def doPlot(self, model=None):
        if model.results is None:
            return
        ax = self.ax
        ax1=self.ax1
        # plot result of odeint solver
        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if ax == None:
            ax = plt.subplot()
            ax1=ax.twinx()
            QTPlotting = False  # actually, we are just using CLI and showing the plot
        ax.clear()
        ax1.clear()
        t=model.t
        ycar = model.results[:,0]
        ywheel=model.results[:,2]
        accel=model.accel

        if self.chk_LogX.isChecked():
            ax.set_xlim(0.001,model.tmax)
            ax.set_xscale('log')
        else:
            ax.set_xlim(0.0, model.tmax)
            ax.set_xscale('linear')

        if self.chk_LogY.isChecked():
            ax.set_ylim(0.0001,max(ycar.max(), ywheel.max()*1.05))
            ax.set_yscale('log')
        else:
            ax.set_ylim(0.0, max(ycar.max(), ywheel.max()*1.05))
            ax.set_yscale('linear')

        ax.plot(t, ycar, 'b-', label='Body Position')
        ax.plot(t, ywheel, 'r-', label='Wheel Position')
        if self.chk_ShowAccel.isChecked():
            ax1.plot(t, accel, 'g-', label='Body Accel')
            ax1.axhline(y=accel.max(), color='orange')  # horizontal line at accel.max()
            ax1.set_yscale('log' if self.chk_LogAccel.isChecked() else 'linear')

        # add axis labels
        ax.set_ylabel("Vertical Position (m)", fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel("time (s)", fontsize='large' if QTPlotting else 'medium')
        ax1.set_ylabel("Y'' (g)", fontsize = 'large' if QTPlotting else 'medium')
        ax.legend()

        ax.axvline(x=model.tramp)  # vertical line at tramp
        ax.axhline(y=model.ymag)  # horizontal line at ymag
        # modify the tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks
        ax1.tick_params(axis='both', which='both', direction='in', right=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks
        # show the plot
        if QTPlotting == False:
            plt.show()
        else:
            self.canvas.draw()


class CarController():
    def __init__(self, args):
        """
        This is the controller I am using for the quarter car model.
        """
        self.input_widgets, self.display_widgets = args
        #unpack widgets with same names as they have on the GUI
        self.le_m1, self.le_v, self.le_k1, self.le_c1, self.le_m2, self.le_k2, self.le_ang, \
         self.le_tmax, self.chk_IncludeAccel = self.input_widgets

        self.gv_Schematic, self.chk_LogX, self.chk_LogY, self.chk_LogAccel, \
        self.chk_ShowAccel, self.lbl_MaxMinInfo, self.layout_horizontal_main = self.display_widgets

        self.model = CarModel()
        self.view = CarView(args)

        #self.chk_IncludeAccel=qtw.QCheckBox()

    def ode_system(self, X, t):
        # define the forcing function equation for the linear ramp
        # It takes self.tramp time to climb the ramp, so y position is
        # a linear function of time.
        if t < self.model.tramp:
            y = self.model.ymag * (t / self.model.tramp)
        else:
            y = self.model.ymag

        x1 = X[0]  # car position in vertical direction
        x1dot = X[1]  # car velocity  in vertical direction
        x2 = X[2]  # wheel position in vertical direction
        x2dot = X[3]  # wheel velocity in vertical direction

        # write the non-trivial equations in vertical direction
        x1ddot = (self.model.k2 * (x2 - x1) - self.model.c1 * (x2dot - x1dot) + self.model.k1 * (0 - x1)) / self.model.m1
        x2ddot = (self.model.k2 * (y - x2) - self.model.k2 * (x2 - x1) + self.model.c1 * (x2dot - x1dot)) / self.model.m2

        # return the derivatives of the input state vector
        return [x1dot, x1ddot, x2dot, x2ddot]

    def calculate(self, doCalc=True):
        """
        I will first set the basic properties of the car model and then calculate the result
        in another function doCalc.
        """
        #Step 1.  Read from the widgets
        self.model.m1 = float(self.le_m1.text())
        self.model.m2 = float(self.le_m2.text())
        self.model.c1 = float(self.le_c1.text())
        self.model.k1 = float(self.le_k1.text())
        self.model.k2 = float(self.le_k2.text())
        self.model.v = float(self.le_v.text())

        # Recalculate min and max k values based on some predefined method or assumptions
        self.mink1 = 28947.64  # minimum spring constant for suspension in N/m
        self.maxk1 = 57895.28  # maximum spring constant for suspension in N/m
        self.mink2 = 115790.55  # minimum spring constant for tire in N/m
        self.maxk2 = 231581.1  # maximum spring constant for tire in N/m

        ymag=6.0/(12.0*3.3)   #This is the height of the ramp in m
        if ymag is not None:
            self.model.ymag = ymag
        self.model.yangdeg = float(self.le_ang.text())
        self.model.tmax = float(self.le_tmax.text())
        if(doCalc):
            self.doCalc()
        self.SSE((self.model.k1, self.model.c1, self.model.k2), optimizing=False)
        self.view.updateView(self.model)

    def setWidgets(self, w):
        self.view.setWidgets(w)
        self.chk_IncludeAccel=self.view.chk_IncludeAccel

    def doCalc(self, doPlot=True, doAccel=True):
        """
        This solves the differential equations for the quarter car model.
        :param doPlot:
        :param doAccel:
        :return:
        """
        v = 1000 * self.model.v / 3600  # convert speed to m/s from kph
        self.model.angrad = self.model.yangdeg * math.pi / 180.0  # convert angle to radians
        self.model.tramp = self.model.ymag / (math.sin(self.model.angrad) * v)  # calculate time to traverse ramp

        self.model.t=np.linspace(0,self.model.tmax,2000)
        ic = [0, 0, 0, 0]
        # run odeint solver
        self.model.results = odeint(self.ode_system, ic, self.model.t)
        if doAccel:
            self.calcAccel()
        if doPlot:
            self.doPlot()

    def calcAccel(self):
        """
        Calculate the acceleration in the vertical direction using the forward difference formula.
        """
        N=len(self.model.t)
        self.model.accel=np.zeros(shape=N)
        vel=self.model.results[:,1]
        for i in range(N):
            if i==N-1:
                h = self.model.t[i] - self.model.t[i-1]
                self.model.accel[i]=(vel[i]-vel[i-1])/(9.81*h)  # backward difference of velocity
            else:
                h = self.model.t[i + 1] - self.model.t[i]
                self.model.accel[i] = (vel[i + 1] - vel[i]) / (9.81 * h)  # forward difference of velocity
            # else:
            #     self.model.accel[i]=(vel[i+1]-vel[i-1])/(9.81*2.0*h)  # central difference of velocity
        self.model.accelMax=self.model.accel.max()
        return True

    def OptimizeSuspension(self):
        """
        Step 1:  set parameters based on GUI inputs by calling self.set(doCalc=False)
        Step 2:  make an initial guess for k1, c1, k2
        Step 3:  optimize the suspension
        :return:
        """
        #Step 1:

        self.calculate(doCalc=False)
        #Step 2:

        x0= np.array([self.model.k1, self.model.c1, self.model.k2])
        #Step 3:

        bounds = [(self.model.mink1, self.model.maxk1), (100, 2000), (self.model.mink2, self.model.maxk2)]

        print("Initial guesses:", x0)
        print("Bounds:", bounds)

        answer= minimize(self.SSE, x0, method='Nelder-Mead', bounds=bounds) #use the Nelder-Mead method to minimize the SSE function (our objective function)

        self.model.k1, self.model.c1, self.model.k2 = answer.x  # Update model with optimized values
        self.view.updateView(self.model)  # Refresh the GUI to show new results

    def SSE(self, vals, optimizing=True):
        """
        Calculates the sum of square errors between the contour of the road and the car body.
        :param vals:
        :param optimizing:
        :return:
        """
        k1, c1, k2=vals  #unpack the new values for k1, c1, k2
        self.model.k1=k1
        self.model.c1=c1
        self.model.k2=k2
        self.doCalc(doPlot=False)  #solve the odesystem with the new values of k1, c1, k2
        SSE=0
        for i in range(len(self.model.results[:,0])):
            t=self.model.t[i]
            y=self.model.results[:,0][i]
            if t < self.model.tramp:
                ytarget = self.model.ymag * (t / self.model.tramp)
            else:
                ytarget = self.model.ymag
            SSE+=(y-ytarget)**2

        #some penalty functions if the constants are too small
        if optimizing:
            if k1<self.model.mink1 or k1>self.model.maxk1:
                SSE+=100
            if c1<10:
                SSE+=100
            if k2<self.model.mink2 or k2>self.model.maxk2:
                SSE+=100

            # I'm overlaying a gradient in the acceleration limit that scales with distance from a target squared.
            if self.model.accelMax>self.model.accelLim and self.chk_IncludeAccel.isChecked():
                # need to soften suspension
                SSE+=(self.model.accelMax-self.model.accelLim)**2
        self.model.SSE=SSE
        return SSE

    def doPlot(self):
        self.view.doPlot(self.model)
#endregion
#endregion

def main():
    QCM = CarController()
    QCM.doCalc()

if __name__ == '__main__':
    main()
