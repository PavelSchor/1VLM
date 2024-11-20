import pygame
import numpy as np

class JoystickController(object):
	def __init__(self):
		

		self.iAileron=0
		self.iElevator=1
		self.iRudder=2
		self.iThrottle=3

		self.dAileron=0.0
		self.dElevator=0.0
		self.dRudder=0.0
		self.dThrottle=0.0
		self.inputs=[0.0,0.0,0.0,0.0]
		
		self.maxAileron=1.0
		self.maxElevator=1.0
		self.maxRudder=1.0
		self.maxThrottle=1.0
		
		self.minimum=np.zeros(4)
		self.maximum=np.zeros(4)
		
		
	def initJoysticks(self):
		pygame.init()
		self.clock = pygame.time.Clock()
		pygame.joystick.init()
		joystick_count = pygame.joystick.get_count()

    
		 #for each joystick:
			#for i in range(joystick_count):
		self.joystick = pygame.joystick.Joystick(0)
		self.joystick.init()
	
	def calibrate(self,t=5.):
		axes = self.joystick.get_numaxes()
		pygame.event.pump()
		
		start_ticks=pygame.time.get_ticks() #starter tick
		while True: # mainloop
			seconds=(pygame.time.get_ticks()-start_ticks)/1000 #calculate how many seconds
			if seconds>t: # if more than 10 seconds close the game
				break
			else:
				pygame.event.pump()
				for i in range( axes ):
					self.inputs[i] = self.joystick.get_axis( i )
					if self.inputs[i]>self.maximum[i]:
						self.maximum[i]=self.inputs[i]
					if self.inputs[i]<self.minimum[i]:
						self.minimum[i]=self.inputs[i]
				pygame.time.delay(100)
	def capture(self):
		axes = self.joystick.get_numaxes()
		pygame.event.pump()
		for i in range( axes ):
			self.inputs[i] = self.joystick.get_axis( i )
			if self.inputs[i]>0.0:
				self.inputs[i]/=self.maximum[i]*np.sign(self.inputs[i])
			if self.inputs[i]<0.0:
				self.inputs[i]/=self.minimum[i]*np.sign(self.inputs[i])
			
			#self.inputs[i] = (self.joystick.get_axis( i ) - self.minimum[i])/(self.maximum[i] - self.minimum[i])
			#textPrint.pprint(screen, "Axis {} value: {:>6.3f}".format(i, axis) )
			pass
		#self.inputs[0]=0.0;self.inputs[2]=0.0;self.inputs[1]=1.0;
	
	def applyInputs(self,result):
		self.capture()
		result[self.maskAileronL]= self.maxAileron*self.inputs[self.iAileron]
		result[self.maskAileronR]=-self.maxAileron*self.inputs[self.iAileron]

		result[self.maskElevator]=self.maxElevator*self.inputs[self.iElevator]*-1
		result[self.maskRudder]=self.maxRudder*self.inputs[self.iRudder]*-1
		
## -------- Main Program Loop ---
#--------
#while done==False:
    ## EVENT PROCESSING STEP
    #for event in pygame.event.get(): # User did something
        #if event.type == pygame.QUIT: # If user clicked close
            #done=True # Flag that we are done so we exit this loop
        
        ## Possible joystick actions: JOYAXISMOTION JOYBALLMOTION JOYBUTTONDOWN JOYBUTTONUP JOYHATMOTION
        #if event.type == pygame.JOYBUTTONDOWN:
            #print("Joystick button pressed.")
        #if event.type == pygame.JOYBUTTONUP:
            #print("Joystick button released.")
            
 
    ## DRAWING STEP
    ## First, clear the screen to white. Don't put other drawing commands
    ## above this, or they will be erased with this command.
    #screen.fill(WHITE)
    #textPrint.reset()

    ## Get count of joysticks

    
        #textPrint.pprint(screen, "Joystick {}".format(i) )
        #textPrint.indent()
    
        ## Get the name from the OS for the controller/joystick
        #name = joystick.get_name()
        #textPrint.pprint(screen, "Joystick name: {}".format(name) )
        
        ## Usually axis run in pairs, up/down for one, and left/right for
        ## the other.
        #axes = joystick.get_numaxes()
        #textPrint.pprint(screen, "Number of axes: {}".format(axes) )
        #textPrint.indent()
        
        #for i in range( axes ):
            #axis = joystick.get_axis( i )
            #textPrint.pprint(screen, "Axis {} value: {:>6.3f}".format(i, axis) )
        #textPrint.unindent()
            
        #buttons = joystick.get_numbuttons()
        #textPrint.pprint(screen, "Number of buttons: {}".format(buttons) )
        #textPrint.indent()

        #for i in range( buttons ):
            #button = joystick.get_button( i )
            #textPrint.pprint(screen, "Button {:>2} value: {}".format(i,button) )
        #textPrint.unindent()
            
        ## Hat switch. All or nothing for direction, not like joysticks.
        ## Value comes back in an array.
        #hats = joystick.get_numhats()
        #textPrint.pprint(screen, "Number of hats: {}".format(hats) )
        #textPrint.indent()

        #for i in range( hats ):
            #hat = joystick.get_hat( i )
            #textPrint.pprint(screen, "Hat {} value: {}".format(i, str(hat)) )
        #textPrint.unindent()
        
        #textPrint.unindent()

    
    ## ALL CODE TO DRAW SHOULD GO ABOVE THIS COMMENT
    
    ## Go ahead and update the screen with what we've drawn.
    #pygame.display.flip()

    ## Limit to 20 frames per second
    #clock.tick(20)
    
## Close the window and quit.
## If you forget this line, the program will 'hang'
## on exit if running from IDLE.
#pygame.quit ()