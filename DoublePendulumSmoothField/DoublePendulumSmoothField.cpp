#define NOMINMAX 1
#include <Windows.h>
#include <algorithm>

#include <optional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include <thread>
#include <locale>

#include "SAFGUIF/SAFGUIF.h"
#include "multidimentional_point.h"

constexpr int list_size = 1024*32;

struct pendulum {
	using pt4 = dixelu::point<4>;

	double g;
	double length, mass;
	double theta1, theta2, 
		p1, p2;
	double x, y;
	pt4 color;

	pt4 get_drivative_at(pt4 point) const {
		double o_theta1 = point[0];
		double o_theta2 = point[1];
		double ptheta1 = point[2];
		double ptheta2 = point[3];

		double theta1dot = (6. / (mass * length * length)) * (2. * ptheta1 - 3. * ptheta2 * std::cos(o_theta1 - o_theta2)) / (16. - 9. * std::pow(std::cos(o_theta1 - o_theta2), 2.));
		double theta2dot = (6. / (mass * length * length)) * (8. * ptheta2 - 3. * ptheta1 * std::cos(o_theta1 - o_theta2)) / (16. - 9. * std::pow(std::cos(o_theta1 - o_theta2), 2.));

 		double p1dot = -0.5 * (mass * length * length) * (theta1dot * theta2dot * std::sin(o_theta1 - o_theta2) + 3 * (g / length) * std::sin(o_theta1));
		double p2dot = -0.5 * (mass * length * length) * (-1. * theta1dot * theta2dot * std::sin(o_theta1 - o_theta2) + (g / length) * std::sin(o_theta2));

		return pt4({ theta1dot, theta2dot, p1dot, p2dot });
	}
	inline pt4 get_params() const {
		return pt4{ theta1, theta2, p1, p2 };
	}
	inline void set_params(pt4 point) {
		theta1 = point[0];
		theta2 = point[1];
		p1 = point[2];
		p2 = point[3];
	}

	inline void evaluate(double h = 0.01) {
		auto params = get_params();

		auto k1 = get_drivative_at(params);
		auto k2 = get_drivative_at(params + 0.5 * h * k1);
		auto k3 = get_drivative_at(params + 0.5 * h * k2);
		auto k4 = get_drivative_at(params + h * k3);

		params += (h / 6) * (k1 + 2. * k2 + 2. * k3 + k4);
		set_params(params);
	}

	void draw() {
		double x1 = x + length * std::sin(theta1);
		double x2 = x1 + length * std::sin(theta2);
		double y1 = y - length * std::cos(theta1);
		double y2 = y1 - length * std::cos(theta2);

		glColor4d(color[0], color[1], color[2], color[3]);
		glBegin(GL_LINE_STRIP);
		glVertex2d(x, y);
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glEnd();
	}
};

auto hsva_to_rgba = [](pendulum::pt4 hsva)->pendulum::pt4 {
	double      hh, p, q, t, ff;
	int i;

	if (hsva[1] <= 0.0)        // < is bogus, just shuts up warnings
		return { hsva[2], hsva[2], hsva[2], hsva[4] };

	hh = hsva[0];
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = hsva[2] * (1.0 - hsva[1]);
	q = hsva[2] * (1.0 - (hsva[1] * ff));
	t = hsva[2] * (1.0 - (hsva[1] * (1.0 - ff)));

	switch (i) {
	case 0:
		return { hsva[2], t, p, hsva[4] };
	case 1:
		return { q, hsva[2], p, hsva[4] };
	case 2:
		return { p, hsva[2], t, hsva[4] };
	case 3:
		return { p, q, hsva[2], hsva[4] };
	case 4:
		return { t, p, hsva[2], hsva[4] };
	default:
		break;
	}
	return { hsva[2], p, q, hsva[4] };
};

void check_for_interpolation(
	std::list<pendulum>& pendulums,
	size_t cur_index,
	std::list<pendulum>::iterator cur_it
	) {
	if (pendulums.size() < 6)
		throw std::runtime_error("Less than 6 pendulums");

	constexpr double angle_eps = 0.01;
	auto step_it = cur_it;

	if (cur_index >= 1) {
		step_it--;
		if (std::abs(step_it->theta1 - cur_it->theta1) < angle_eps &&
			std::abs(step_it->theta2 - cur_it->theta2) < angle_eps)
			return;
		auto old_step = step_it;
		step_it = cur_it;
		cur_it = old_step;
	}
	else {
		return;
	}

	if (!cur_index || cur_index >= pendulums.size() - 2)
		return;

	auto new_it = cur_it;
	new_it--;
	auto p0 = (new_it)->get_params();
	new_it++;
	auto p1 = (new_it)->get_params();
	new_it++;
	auto p2 = (new_it)->get_params();
	new_it++;
	auto p3 = (new_it)->get_params();


	double max_diff = std::max(std::abs(p1[0] - p2[0]), std::abs(p1[1] - p2[1]));
	auto n_prob = max_diff / angle_eps;
	int n = (std::floor(n_prob));
	if (RANDPOSFLOAT(1) < (n_prob - n)) 
		n++;

	auto interpolate = [&p0, &p1, &p2, &p3](double x) -> pendulum::pt4 {
		x *= angle_eps;
		return
			p0 * (x * (x - angle_eps) * (x - 2 * angle_eps) / (-6 * std::pow(angle_eps, 3))) +
			p1 * ((x * x - angle_eps * angle_eps) * (x - 2 * angle_eps) / (2 * std::pow(angle_eps, 3)))  +
			p2 * (x * (x + angle_eps) * (x - 2 * angle_eps) / (-2 * std::pow(angle_eps, 3))) +
			p3 * (x * (x * x - angle_eps * angle_eps) / (6 * std::pow(angle_eps, 3)));
	};

	auto simple_interpolation = [&p1, &p2](double x) -> pendulum::pt4 {
		return p1 * x + p2 * (1 - x);
	};

	for (int i = 1; i < n; i++) {
		double alpha =  double(i) / n;
		pendulum pd = *cur_it;
		pd.set_params(simple_interpolation(alpha));
		//pd.color = cur_it->color * alpha + step_it->color * (1 - alpha);
		pendulums.insert(step_it, pd);
	}

	if (cur_index > pendulums.size() / 2) 
		while (pendulums.size() > list_size) 
			pendulums.pop_front();
	else 
		while (pendulums.size() > list_size) 
			pendulums.pop_back();
}

void evalute_set_of_pendulums(std::list<pendulum>& pendulums) {
	size_t pend_index = 0;
	for (auto it = pendulums.begin(); it != pendulums.end(); it++, pend_index++) {
		it->evaluate(0.0125);
		check_for_interpolation(pendulums, pend_index, it);
	}
}


void Init() {
	static ButtonSettings* BS_List_Black_Small = new ButtonSettings(System_White, 0, 0, 100, 5, 1, 0, 0, 0xFFEFDFFF, 0x00003F7F, 0x7F7F7FFF);

	MoveableWindow* T = new MoveableResizeableWindow("Main window", System_White, -200, 200, 400, 400, 0xFF, 0x3F3F3F7F, 0);
	((MoveableResizeableWindow*)T)->AssignMinDimentions(200, 200);


	(*WH)["MAIN"] = T;

	WH->EnableWindow("MAIN");

	//__main();
}

///////////////////////////////////////
/////////////END OF USE////////////////
///////////////////////////////////////

#define _WH(MainWindow,Element) ((*(*WH)[MainWindow])[Element])//...uh
#define _WH_t(MainWindow,Element,Type) ((Type)_WH(MainWindow,Element))

void onTimer(int v);
void mDisplay() {
	static std::list<pendulum> pends;

	glClear(GL_COLOR_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	auto color_by_index = [](double i) -> decltype(pendulum::color) {
		return { 0.125*double(i) / list_size, 0.5 + double(i) * 0.5 / list_size, 1., 0.05 };
	};

	if (FIRSTBOOT) {
		FIRSTBOOT = 0;

		for (int i = 0; i < list_size; i++) {
			pendulum pd;
			const double w = 100;
			pd.x = 0;//(0.5 - double(i)/list_size)*2.*w;
			pd.y = 0;
			pd.color = color_by_index(i);
			pd.g = 10;
			pd.length = 75;
			pd.mass = 10;
			pd.p1 = 0;
			pd.p2 = 0;
			pd.theta1 = 2.9;
			pd.theta2 = 1.3 + (1e-3 / list_size) * (double(i) + 1);
			pends.push_back(pd);
		}

		//ANIMATION_IS_ACTIVE = !ANIMATION_IS_ACTIVE;
		onTimer(0);
	}

	if(ANIMATION_IS_ACTIVE)
		evalute_set_of_pendulums(pends);

	int cnt = 0;
	for (auto& p : pends) {
		cnt++;
		p.color = color_by_index(cnt);
		p.draw();
	}

	glutSwapBuffers();
}

void mInit() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D((0 - RANGE) * (WindX / WINDXSIZE), RANGE * (WindX / WINDXSIZE), (0 - RANGE) * (WindY / WINDYSIZE), RANGE * (WindY / WINDYSIZE));
}

void onTimer(int v) {
	glutTimerFunc(16, onTimer, 0);
	if (true || ANIMATION_IS_ACTIVE) {
		mDisplay();
		++TimerV;
	}
}

void OnResize(int x, int y) {
	WindX = x;
	WindY = y;
	mInit();
	glViewport(0, 0, x, y);
}
void inline rotate(float& x, float& y) {
	float t = x * cos(ROT_RAD) + y * sin(ROT_RAD);
	y = 0.f - x * sin(ROT_RAD) + y * cos(ROT_RAD);
	x = t;
}
void inline absoluteToActualCoords(int ix, int iy, float& x, float& y) {
	float wx = WindX, wy = WindY;
	x = ((float)(ix - wx * 0.5f)) / (0.5f * (wx / (RANGE * (WindX / WINDXSIZE))));
	y = ((float)(0 - iy + wy * 0.5f)) / (0.5f * (wy / (RANGE * (WindY / WINDYSIZE))));
	rotate(x, y);
}
void mMotion(int ix, int iy) {
	float fx, fy;
	absoluteToActualCoords(ix, iy, fx, fy);
	MXPOS = fx;
	MYPOS = fy;
	if (WH)WH->MouseHandler(fx, fy, 0, 0);
}
void mKey(BYTE k, int x, int y) {
	if (WH)WH->KeyboardHandler(k);

	if (k == 's') 
		ANIMATION_IS_ACTIVE = !ANIMATION_IS_ACTIVE;

	if (k == 27)
		exit(1);
}
void mClick(int butt, int state, int x, int y) {
	float fx, fy;
	CHAR Button, State = state;
	absoluteToActualCoords(x, y, fx, fy);
	Button = butt - 1;
	if (state == GLUT_DOWN)State = -1;
	else if (state == GLUT_UP)State = 1;
	if (WH)WH->MouseHandler(fx, fy, Button, State);
}
void mDrag(int x, int y) {
	mMotion(x, y);
}
void mSpecialKey(int Key, int x, int y) {
	auto modif = glutGetModifiers();
	if (!(modif & GLUT_ACTIVE_ALT)) {
		switch (Key) {
		case GLUT_KEY_DOWN:		if (WH)WH->KeyboardHandler(1);
			break;
		case GLUT_KEY_UP:		if (WH)WH->KeyboardHandler(2);
			break;
		case GLUT_KEY_LEFT:		if (WH)WH->KeyboardHandler(3);
			break;
		case GLUT_KEY_RIGHT:	if (WH)WH->KeyboardHandler(4);
			break;
		}
	}
	if (modif == GLUT_ACTIVE_ALT && Key == GLUT_KEY_DOWN) {
		RANGE *= 1.1;
		OnResize(WindX, WindY);
	}
	else if (modif == GLUT_ACTIVE_ALT && Key == GLUT_KEY_UP) {
		RANGE /= 1.1;
		OnResize(WindX, WindY);
	}
}
void mExit(int a) {

}

int main(int argc, char** argv) {
	std::ios_base::sync_with_stdio(false);//why not
#ifdef _DEBUG 
	ShowWindow(GetConsoleWindow(), SW_SHOW);
#else // _DEBUG 
	ShowWindow(GetConsoleWindow(), SW_HIDE);
#endif
	ShowWindow(GetConsoleWindow(), SW_SHOW);

	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
	//srand(1);
	//srand(clock());
	InitASCIIMap();
	//cout << to_string((WORD)0) << endl;

	srand(TIMESEED());
	__glutInitWithExit(&argc, argv, mExit);
	//cout << argv[0] << endl;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_MULTISAMPLE);
	glutInitWindowSize(WINDXSIZE, WINDYSIZE);
	//glutInitWindowPosition(50, 0);
	glutCreateWindow(WINDOWTITLE);

	hWnd = FindWindowA(NULL, WINDOWTITLE);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE);//_MINUS_SRC_ALPHA
	glEnable(GL_BLEND);

	//glEnable(GL_POLYGON_SMOOTH);//laggy af
	glEnable(GL_LINE_SMOOTH);//GL_POLYGON_SMOOTH
	glEnable(GL_POINT_SMOOTH);

	glShadeModel(GL_SMOOTH);
	//glEnable(GLUT_MULTISAMPLE);

	glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);//GL_FASTEST//GL_NICEST
	glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

	glutMouseFunc(mClick);
	glutReshapeFunc(OnResize);
	glutSpecialFunc(mSpecialKey);
	glutMotionFunc(mDrag);
	glutPassiveMotionFunc(mMotion);
	glutKeyboardFunc(mKey);
	glutDisplayFunc(mDisplay);
	mInit();
	glutMainLoop();
	return 0;
}