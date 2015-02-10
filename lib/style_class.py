
import matplotlib

class style_container():
    def __init__(self, istyle = 'Plain', useRoot = False, addplots = ['', '', ''], addheights = [0, 0, 0], cmsPositon = "upper right", legendPosition = "upper right"):
        self._style = istyle

        self._useRoot = useRoot

        self.addplots = addplots
        self.addheights = addheights

        if self._useRoot:
            self.Set_Root_style()
        else:
            matplotlib.rcParams.update({'font.size': 10})
            matplotlib.rcParams.update({'lines.linewidth' : 1})
        #rc('text', usetex=True)
        # self._xaxis_title      = self._allHists[0].xaxis.GetTitle()
        # self._yaxis_title      = self._allHists[0].yaxis.GetTitle()
        self._xaxis_title      = 'bla'
        self._yaxis_title      = '#epsilon'
        self._additional_text  = 'Preliminary'
        self._y_label_offset   = -0.11
        self._error_bands_ecol = ['darkmagenta','darkcyan']
        self._error_bands_fcol = ['m','cyan']
        self._error_bands_alph = 0.7
        self._error_bands_labl = ['Sys. uncert. 1','Sys. uncert. 2']
        self._error_bands_center = 'ref'
        self._error_stacking = 'No'
        self._spine_line_width = 0.5
        self._logx = False
        self._logy = True
        self._ymin = -1
        self._ymax = -1
        self._xmin = -1
        self._xmax = -1
        self.cmsTextPosition=position(cmsPositon,isText=True)
        self.LegendPosition=position(legendPosition)
        if self._style == 'CMS':
            self.Set_CMS_style()
        elif self._style == 'Plain':
            self.Set_Plain_style()
        elif self._style == 'Cool':
            self.Set_Cool_style()

    def __del__(self):
        pass

    def Set_CMS_style(self):
        self._add_cms_text           = True
        self._add_lumi_text          = True
        self._label_text_color       = 'black'
        self._annotation_text_color  = 'black'
        self._bg_color               = 'w'
        self._ref_line_color         = 'blue'
        self._spine_color            = 'black'
        self._tick_color             = 'black'
        self._marker_style           = 'o'
        self._marker_size            = 3
        self._marker_color           = 'black'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'row'
        self._show_minor_tick_labels = False
        self._legend_font_size       = 9
        if self.addplots[0] != '':
            self.cmsTextPosition.addYspace(  -0.9 * self.addheights[0] / 100.)
        if self.addplots[1] != '':
            self.cmsTextPosition.addYspace(  0.9 * self.addheights[1] / 100.)
        if self.addplots[2] != '':
            self.cmsTextPosition.addYspace(  0.9 * self.addheights[2] / 100.)

    def Set_Plain_style(self):
        self._add_cms_text           = False
        self._add_lumi_text          = False
        self._label_text_color       = 'black'
        self._annotation_text_color  = 'black'
        self._bg_color               = 'w'
        self._ref_line_color         = 'blue'
        self._spine_color            = 'black'
        self._tick_color             = 'black'
        self._marker_style           = 'o'
        self._marker_size            = 4
        self._marker_color           = 'black'
        self._marker_error_cap_width = 1
        self._cms_text_alignment     = 'row'
        self._show_minor_tick_labels = True
        self._legend_font_size       = 10
        if self.addplots[0] != '':
            self.cmsTextPosition.addYspace(  -0.8 * self.addheights[0] / 100.)
        if self.addplots[1] != '':
            self.cmsTextPosition.addYspace(  0.8 * self.addheights[1] / 100.)
        if self.addplots[2] != '':
            self.cmsTextPosition.addYspace(  0.8 * self.addheights[2] / 100.)

    def Set_Cool_style(self):
        self._add_cms_text           = True
        self._add_lumi_text          = True
        self._label_text_color       = 'white'
        self._annotation_text_color  = 'white'
        self._bg_color               = '#07000d'
        self._ref_line_color         = 'y'
        self._spine_color            = '#5998ff'
        self._tick_color             = 'w'
        self._marker_style           = 'o'
        self._marker_size            = 3
        self._marker_color           = 'lightgray'
        self._marker_error_cap_width = 0
        self._cms_text_alignment     = 'column'
        self._show_minor_tick_labels = False
        self._legend_font_size       = 9

    def Set_Root_style(self):
        self.cmsTextFont          = 61   # Fonts
        self.lumiTextFont         = 42
        self.extraTextFont        = 52
        self.additionalTextFont   = 42
        self.cmsTextSize          = 0.9  #Text sizes
        self.lumiTextSize         = 0.6
        self.extraTextSize        = 0.76*self.cmsTextSize
        self.additionalTextSize   = 1.0*self.extraTextSize
        self.legendTextSize       = self.extraTextSize*0.8
        self.lumiTextOffset       = 0.2
        self.extraTextOffset      = 2.5  # only used in outOfFrame version
        self.axisTextSize         = 0.9
        self.axisOffset           = 1.3
        self._ratio_pad           ={}
        self.rootMemory           =[]

    ## Function to set properties of the plotting axis
    #
    # This function sets axis properties like the y-range or
    # if any axis should be logarithmic.
    # @param[in] logx Boolean if the x-axis should be logarithmic (Default = False)
    # @param[in] logy Boolean if the y-axis should be logarithmic (Default = True)
    # @param[in] ymin Minimum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] ymax Maximum plotting range for the y-axis (Default = -1 automatic values)
    # @param[in] xmin Minimum plotting range for the x-axis (Default = -1 range from hist)
    # @param[in] xmax Maximum plotting range for the x-axis (Default = -1 range from hist)
    def Set_axis(self, logx = False, logy = True, ymin = -1, ymax = -1, xmin = -1, xmax = -1):
        self._logx = logx
        self._logy = logy
        self._ymin = ymin
        self._ymax = ymax
        self._xmin = xmin
        self._xmax = xmax

class position():
    def __init__(self,positiontext="upper right", refference="", isText=False):

        self._positiontext=positiontext
        if not isinstance(positiontext,str):
            self._definedCoorinates=True
            self._valign="left"
            self._align="left"
        else:
            self._definedCoorinates=False
            self._valign=self._positiontext.split(" ")[0]
            self._align=self._positiontext.split(" ")[1]
        self.addY=0
        self.addX=0
        self._isText=isText
        self._correctcms={"left":0.,
                    "middle":0.,
                    "right":-0.15,
                    "upper":-0.04,
                    "center":0.,
                    "lower":0.,
        }
        if self._definedCoorinates:
            self._x=self._positiontext[0]
            self._y=self._positiontext[1]

    def __eq__(self,other):
        return (self._positiontext==other._positiontext)

    def addYspace(self,y):
        if self._valign=="upper" and y<0.:
            self.addY+=y
        elif self._valign!="upper" and self.getY()+y>0.1:
            self.addY+=y

    def addXspace(self,x):
        self.addX+=x



    def setPosition(self,positiontext):
        self._positiontext=positiontext
        self.valign=self._positiontext.split(" ")[0]
        self.align=self._positiontext.split(" ")[1]

    def getText(self):
        return self._positiontext

    def getX(self):
        if self._definedCoorinates:
            return self._x
        alignDict={
                    "left":0.12,
                    "middle":0.5,
                    "right":0.95,
        }
        if self._isText:
            return self.addX+alignDict[self._align]+self._correctcms[self._align]
        return self.addX+alignDict[self._align]

    def getY(self):
        if self._definedCoorinates:
            return self._y
        alignDict={
                    "upper":0.95,
                    "center":0.5,
                    "lower":0.12,
        }
        if self._isText:
            return self.addY+alignDict[self._valign]+self._correctcms[self._valign]
        return self.addY+alignDict[self._valign]
