#pragma once
#include "Hamiltionian.h"

std::vector<std::vector<double>> eigenVec_num, eigenVec_analitycal;
std::vector<double> eigenval_num, eigenVal_analityc, Quantum_Well;

namespace ParabolicQW {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::DataVisualization::Charting::Chart^ chart1;
	protected:
	private: System::Windows::Forms::Button^ button1;
	private: System::Windows::Forms::NumericUpDown^ numericUpDown_x_0;
	private: System::Windows::Forms::NumericUpDown^ numericUpDown_width;


	private: System::Windows::Forms::NumericUpDown^ numericUpDown_dx;
	private: System::Windows::Forms::NumericUpDown^ numericUpDown_omega;
	private: System::Windows::Forms::NumericUpDown^ numericUpDown_n;






	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::Label^ label2;

	private: System::Windows::Forms::Label^ label4;
	private: System::Windows::Forms::Label^ label5;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::CheckBox^ checkBox_num;
	private: System::Windows::Forms::CheckBox^ checkBox_analityc;
	private: System::Windows::Forms::Button^ button2;


	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataVisualization::Charting::ChartArea^ chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^ legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^ series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->numericUpDown_x_0 = (gcnew System::Windows::Forms::NumericUpDown());
			this->numericUpDown_width = (gcnew System::Windows::Forms::NumericUpDown());
			this->numericUpDown_dx = (gcnew System::Windows::Forms::NumericUpDown());
			this->numericUpDown_omega = (gcnew System::Windows::Forms::NumericUpDown());
			this->numericUpDown_n = (gcnew System::Windows::Forms::NumericUpDown());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->checkBox_num = (gcnew System::Windows::Forms::CheckBox());
			this->checkBox_analityc = (gcnew System::Windows::Forms::CheckBox());
			this->button2 = (gcnew System::Windows::Forms::Button());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_x_0))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_width))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_dx))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_omega))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_n))->BeginInit();
			this->SuspendLayout();
			// 
			// chart1
			// 
			chartArea1->AxisX->Title = L"x [nm]";
			chartArea1->AxisX->TitleFont = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			chartArea1->AxisY->Title = L"Energy [eV]";
			chartArea1->AxisY->TitleFont = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			chartArea1->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea1);
			legend1->Name = L"Legend1";
			this->chart1->Legends->Add(legend1);
			this->chart1->Location = System::Drawing::Point(280, 31);
			this->chart1->Name = L"chart1";
			series1->BorderWidth = 3;
			series1->ChartArea = L"ChartArea1";
			series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Spline;
			series1->Legend = L"Legend1";
			series1->Name = L"Series1";
			this->chart1->Series->Add(series1);
			this->chart1->Size = System::Drawing::Size(1178, 781);
			this->chart1->TabIndex = 0;
			this->chart1->Text = L"chart1";
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->button1->Location = System::Drawing::Point(12, 31);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(262, 48);
			this->button1->TabIndex = 1;
			this->button1->Text = L"Calculate";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// numericUpDown_x_0
			// 
			this->numericUpDown_x_0->Location = System::Drawing::Point(12, 93);
			this->numericUpDown_x_0->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 100, 0, 0, System::Int32::MinValue });
			this->numericUpDown_x_0->Name = L"numericUpDown_x_0";
			this->numericUpDown_x_0->Size = System::Drawing::Size(75, 20);
			this->numericUpDown_x_0->TabIndex = 2;
			// 
			// numericUpDown_width
			// 
			this->numericUpDown_width->DecimalPlaces = 2;
			this->numericUpDown_width->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 131072 });
			this->numericUpDown_width->Location = System::Drawing::Point(12, 119);
			this->numericUpDown_width->Name = L"numericUpDown_width";
			this->numericUpDown_width->Size = System::Drawing::Size(75, 20);
			this->numericUpDown_width->TabIndex = 3;
			this->numericUpDown_width->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 2, 0, 0, 0 });
			// 
			// numericUpDown_dx
			// 
			this->numericUpDown_dx->DecimalPlaces = 3;
			this->numericUpDown_dx->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 196608 });
			this->numericUpDown_dx->Location = System::Drawing::Point(12, 145);
			this->numericUpDown_dx->Maximum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 0 });
			this->numericUpDown_dx->Minimum = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 196608 });
			this->numericUpDown_dx->Name = L"numericUpDown_dx";
			this->numericUpDown_dx->Size = System::Drawing::Size(75, 20);
			this->numericUpDown_dx->TabIndex = 4;
			this->numericUpDown_dx->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 5, 0, 0, 196608 });
			// 
			// numericUpDown_omega
			// 
			this->numericUpDown_omega->DecimalPlaces = 2;
			this->numericUpDown_omega->Increment = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 131072 });
			this->numericUpDown_omega->Location = System::Drawing::Point(12, 171);
			this->numericUpDown_omega->Name = L"numericUpDown_omega";
			this->numericUpDown_omega->Size = System::Drawing::Size(75, 20);
			this->numericUpDown_omega->TabIndex = 6;
			this->numericUpDown_omega->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 1, 0, 0, 0 });
			// 
			// numericUpDown_n
			// 
			this->numericUpDown_n->Location = System::Drawing::Point(12, 272);
			this->numericUpDown_n->Name = L"numericUpDown_n";
			this->numericUpDown_n->Size = System::Drawing::Size(75, 20);
			this->numericUpDown_n->TabIndex = 7;
			this->numericUpDown_n->Value = System::Decimal(gcnew cli::array< System::Int32 >(4) { 5, 0, 0, 0 });
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label1->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->label1->Location = System::Drawing::Point(93, 95);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(181, 17);
			this->label1->TabIndex = 8;
			this->label1->Text = L"Quantum Well position [nm]";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label2->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->label2->Location = System::Drawing::Point(93, 121);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(164, 17);
			this->label2->TabIndex = 9;
			this->label2->Text = L"Quantum Well width [nm]";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label4->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->label4->Location = System::Drawing::Point(93, 147);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(134, 17);
			this->label4->TabIndex = 11;
			this->label4->Text = L"x-grid sampling [nm]";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label5->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->label5->Location = System::Drawing::Point(93, 173);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(54, 17);
			this->label5->TabIndex = 12;
			this->label5->Text = L"ω [pHz]";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label6->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->label6->Location = System::Drawing::Point(9, 252);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(234, 17);
			this->label6->TabIndex = 13;
			this->label6->Text = L"Number of eigenvalues to calculate:";
			// 
			// checkBox_num
			// 
			this->checkBox_num->AutoSize = true;
			this->checkBox_num->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->checkBox_num->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->checkBox_num->Location = System::Drawing::Point(12, 330);
			this->checkBox_num->Name = L"checkBox_num";
			this->checkBox_num->Size = System::Drawing::Size(178, 21);
			this->checkBox_num->TabIndex = 14;
			this->checkBox_num->Text = L"Calculate eigenvectors\?";
			this->checkBox_num->UseVisualStyleBackColor = true;
			// 
			// checkBox_analityc
			// 
			this->checkBox_analityc->AutoSize = true;
			this->checkBox_analityc->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 10, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->checkBox_analityc->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->checkBox_analityc->Location = System::Drawing::Point(12, 353);
			this->checkBox_analityc->Name = L"checkBox_analityc";
			this->checkBox_analityc->Size = System::Drawing::Size(180, 21);
			this->checkBox_analityc->TabIndex = 15;
			this->checkBox_analityc->Text = L"Analitycal eigenvectors\?";
			this->checkBox_analityc->UseVisualStyleBackColor = true;
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(13, 407);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(261, 307);
			this->button2->TabIndex = 16;
			this->button2->Text = L"SAVE DATA";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackColor = System::Drawing::SystemColors::ActiveCaptionText;
			this->ClientSize = System::Drawing::Size(1479, 824);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->checkBox_analityc);
			this->Controls->Add(this->checkBox_num);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->numericUpDown_n);
			this->Controls->Add(this->numericUpDown_omega);
			this->Controls->Add(this->numericUpDown_dx);
			this->Controls->Add(this->numericUpDown_width);
			this->Controls->Add(this->numericUpDown_x_0);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->chart1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_x_0))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_width))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_dx))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_omega))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->numericUpDown_n))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {
		this->chart1->Series->Clear();

		//Initialize something:
			System::Windows::Forms::DataVisualization::Charting::Series^ series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			series2->ChartArea = L"ChartArea1";
			series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Spline;
			series2->Legend = L"Legend1";
			series2->Name = L"QW";
			series2->BorderWidth = 3;
			this->chart1->Series->Add(series2);

			//double omega = 1; //in pHz 10^15
			double m_0 = 9.10938356; //electron mass in [kg] // 10^-31
			double hbar = 6.582119569; //planck constant in [eV*s] //10^-16
			double e_c = 1.60211; // elementary electric charge //10^-19
			int power = 1; // 10^-31 * (10^15)^2 * 10^(-18) / 10^(-19)

			double x_0 = static_cast<double>(numericUpDown_x_0->Value);
			double width = static_cast<double>(numericUpDown_width->Value);
			double dx = static_cast<double>(numericUpDown_dx->Value);
			double omega = static_cast<double>(numericUpDown_omega->Value);
			int num_of_eigenVal = static_cast<int>(numericUpDown_n->Value);
			//x_1 = -1, x_2 = 1, dx = 0.01; // in nm=10^-9m
			double x_1 = x_0 - width / 2.0;
			double x_2 = x_0 + width / 2.0;

			Hamiltonian Hamil(omega, x_0, x_1, x_2, dx);
		//--------------------------------------------------------------------------------------------------------

		// Build Quantum Well potential with same x-grid as Hamil:
			std::vector<double> potential(static_cast<int>((x_2 - x_1) / dx + 1));
			for (int k = 0; k < potential.size(); k++) {
				double x = x_1 + k * dx;
				potential[k] = 1. / 2. * m_0 * omega * omega * power * (x - x_0) * (x - x_0); // in [eV]
				//this->chart1->Series["QW"]->Points->AddXY(x, potential[k]);
			}
			Quantum_Well = potential; // save to global variable

			double x = x_1 - 0.2;
			while (x <= x_2 + 0.2 + dx){
				if(x <= x_1)
					this->chart1->Series["QW"]->Points->AddXY(x, 1. / 2. * m_0 * omega * omega * power * (x_1 - x_0) * (x_1 - x_0));
				else if (x >= x_2)
					this->chart1->Series["QW"]->Points->AddXY(x, 1. / 2. * m_0 * omega * omega * power * (x_2 - x_0) * (x_2 - x_0));
				else
					this->chart1->Series["QW"]->Points->AddXY(x, 1. / 2. * m_0 * omega * omega * power * (x - x_0) * (x - x_0));
				x += dx;
			}

		//--------------------------------------------------------------------------------------------------------

		//Contacting external executable
			TCHAR Npath[MAX_PATH], Path2[MAX_PATH];
			GetCurrentDirectory(MAX_PATH, Npath);
			SetCurrentDirectory(L"./EXECUTABLE/Test_arma/x64/Debug/");
			GetCurrentDirectory(MAX_PATH, Path2);
			std::ofstream parameters("parameters.dat");

			parameters << x_0 << "\t" << x_1 << "\t" << x_2 << "\t" << dx << "\t" << omega << "\t" << num_of_eigenVal << std::endl;
			parameters.close();
			HINSTANCE hRet = ShellExecute(
				HWND_DESKTOP, //Parent window
				L"open",       //Operation to perform
				L"C:/Users/77swi/source/repos/ParabolicQW/ParabolicQW/EXECUTABLE/Test_arma/x64/Debug/Test_arma.exe",       //Path to program/
				NULL,         //Parameters // -- need to check how to pass parameters
				NULL,         //Default directory
				SW_SHOW);     //How to open

			if ((LONG)hRet <= 32)
			{
				MessageBox::Show("Unable to start program", "Error");
			}
			Sleep(500);
		//--------------------------------------------------------------------------------------------------------


		Hamil.Build_Hamiltonian(potential); // Build Hamiltonian
		int N = potential.size();

		//Check if nr of eigenvalues too big
		if (num_of_eigenVal - 1 >= static_cast<int>( min(potential[0], potential[N-1]) / (hbar * omega * 0.1) - 0.25)) {
			num_of_eigenVal = static_cast<int>(min(potential[0], potential[N - 1]) /(hbar * omega * 0.1 ) - 0.25);
			std::string n = std::to_string(num_of_eigenVal - 1);
			MessageBox::Show("Number of eigenvalues exceeds potential well (numerical results)\n\
Eigenvalues adjusted to\n n_max =  "+ gcnew String(n.c_str()) + "", "Error: Number of E_n");
		}
		// Read eigenvalues from file calculated by external executable
		std::ifstream E("eigenVal.dat");
		std::ifstream eigenVec_file("eigenVec.dat");
		eigenVec_num = std::vector<std::vector<double>>(num_of_eigenVal); // save to global variable
		eigenVec_analitycal = std::vector<std::vector<double>>(num_of_eigenVal);
		for (int k = 0; k < num_of_eigenVal; k++) {
			double energy;
			E >> energy;
			Hamil.eigenVal.push_back(energy);
			eigenval_num.push_back(energy); // save to global variable
			eigenVec_num[k] = std::vector<double>(N);
			for (int m = 0; m < N; m++) {
				eigenVec_file >> eigenVec_num[k][m];
			}
		}
		E.close(); eigenVec_file.close();

		//Hamil.Strum_Ehrlich_Diagonalization(num_of_eigenVal); // well does not work
		//Hamil.Sturm_Diagonalization(2 * num_of_eigenVal); // neither does this
		SetCurrentDirectory(Npath);
		std::string name_base = "E_";
		std::string vec_base = "EigenVec_";
		if (num_of_eigenVal > N) num_of_eigenVal = N; // here put message box
		for (int k = 0; k < num_of_eigenVal; k++) {
			// Adding numerical result
				String^ name = gcnew String((name_base + std::to_string(k)).c_str());
				System::Windows::Forms::DataVisualization::Charting::Series^ series3 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
				series3->ChartArea = L"ChartArea1";
				series3->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
				series3->Legend = L"Legend1";
				series3->Name = name;
				series3->BorderWidth = 3;
				this->chart1->Series->Add(series3);

				double x1 = x_0 - sqrt(2 * Hamil.eigenVal[k] / (m_0 * omega * omega * power));
				double x2 = x_0 + sqrt(2 * Hamil.eigenVal[k] / (m_0 * omega * omega * power));
				this->chart1->Series[name]->Points->AddXY(x1, Hamil.eigenVal[k]);
				this->chart1->Series[name]->Points->AddXY(x2, Hamil.eigenVal[k]);

			//Analitycal result - energy
				name = gcnew String((name_base + std::to_string(k) + "_analitycal").c_str());
				System::Windows::Forms::DataVisualization::Charting::Series^ series4 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
				series4->ChartArea = L"ChartArea1";
				series4->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
				series4->Legend = L"Legend1";
				series4->Name = name;
				series4->BorderWidth = 2;
				this->chart1->Series->Add(series4);
				double E = hbar * omega * 0.1 * (k + 0.5);

				eigenVal_analityc.push_back(E); // save to global variable

				x1 = x_0 - sqrt(2 * E / (m_0 * omega * omega * power));
				x2 = x_0 + sqrt(2 * E / (m_0 * omega * omega * power));
				this->chart1->Series[name]->BorderDashStyle = System::Windows::Forms::DataVisualization::Charting::ChartDashStyle::Dash;
				this->chart1->Series[name]->Points->AddXY(x1, E);
				this->chart1->Series[name]->Points->AddXY(x2, E);

			// Calculating eigenvector for given eigenvalue - analitycal & numerical
				std::vector<double> eigenvector(0);
				if (this->checkBox_analityc->Checked) {
					eigenvector = Hamil.Analitycal_eigenvector(k);

					eigenVec_analitycal[k] = eigenvector; // save to global variable

					name = gcnew String((vec_base + std::to_string(k) + "-analitycal").c_str());
					System::Windows::Forms::DataVisualization::Charting::Series^ series5 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
					series5->ChartArea = L"ChartArea1";
					series5->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
					series5->Legend = L"Legend1";
					series5->Name = name;
					series5->BorderWidth = 2.5;
					this->chart1->Series->Add(series5);
					this->chart1->Series[name]->BorderDashStyle = System::Windows::Forms::DataVisualization::Charting::ChartDashStyle::Dot;
					double maximum = 0;
					for (int m = 0; m < eigenvector.size(); m++) {
						double x = x_1 - 0.2 + m * dx;
						this->chart1->Series[name]->Points->AddXY(x, eigenvector[m]);
						if (eigenvector[m] < maximum)
							maximum = eigenvector[m];
					}
				}
				//------------
				if (this->checkBox_num->Checked) {
					//eigenvector = Hamil.Calculate_eigenvector_Numerov(potential, k);
					//eigenvector = Hamil.Calculate_eigenvector(k);
					eigenvector = std::vector<double>(N);
					double maximal = 0;
					for (int p = 0; p < N; p++) {
						eigenvector[p] = eigenVec_num[k][p];
						if (eigenvector[p] > maximal)
							maximal = eigenvector[p];
					}
					//if (maximum * maximal < 0) maximal *= -1;

					for (int p = 0; p < eigenvector.size(); p++) {
						eigenvector[p] /= maximal;
						if (hbar * omega * 0.1 < 1)
							eigenvector[p] *= (0.5 * hbar * omega * 0.1);
						else
							eigenvector[p] /= (0.5 * hbar * omega * 0.1);
						eigenvector[p] += Hamil.eigenVal[k];
					}

					name = gcnew String((vec_base + std::to_string(k) + "-num").c_str());
					System::Windows::Forms::DataVisualization::Charting::Series^ series6 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
					series6->ChartArea = L"ChartArea1";
					series6->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
					series6->Legend = L"Legend1";
					series6->Name = name;
					series6->BorderWidth = 2.5;
					this->chart1->Series->Add(series6);
					for (int m = 0; m < eigenvector.size(); m++) {
						double x = x_1 + m * dx;
						this->chart1->Series[name]->Points->AddXY(x, eigenvector[m]);
					}
				}
				eigenvector.~vector();

		}

		this->chart1->ChartAreas[0]->RecalculateAxesScale(); //autoscale
		this->chart1->ChartAreas[0]->AxisX->Minimum = x_1 - 0.2;
		this->chart1->ChartAreas[0]->AxisX->Maximum = x_2 + 0.2;
		this->chart1->ChartAreas[0]->AxisY->Minimum = 0;
		this->chart1->ChartAreas[0]->AxisY->Maximum = 1.15*potential[potential.size() - 1];
		//Hamil.~Hamiltonian();
		potential.~vector();
	}

	//Save to file
	private: System::Void button2_Click(System::Object^ sender, System::EventArgs^ e) {
		double x_0 = static_cast<double>(numericUpDown_x_0->Value);
		double width = static_cast<double>(numericUpDown_width->Value);
		double dx = static_cast<double>(numericUpDown_dx->Value);
		double omega = static_cast<double>(numericUpDown_omega->Value);
		int num_of_eigenVal = static_cast<int>(numericUpDown_n->Value);
		double x_1 = x_0 - width / 2.0;
		double x_2 = x_0 + width / 2.0;

		int N = static_cast<int>((x_2 - x_1) / dx + 1);

		SaveFileDialog saveFileDialog;
		TCHAR Npath[MAX_PATH];
		GetCurrentDirectory(MAX_PATH, Npath);
		CreateDirectory(L"./DATA/", NULL);
		SetCurrentDirectory(L"./DATA/");

		saveFileDialog.InitialDirectory = L"./DATA/";
		saveFileDialog.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
		saveFileDialog.FilterIndex = 1;
		saveFileDialog.RestoreDirectory = true;

		std::ostringstream append;
		append << std::setprecision(2) << x_0 << ",width = " << width << ",omega=" << omega;
		std::string name = "QW_Harmonic_oscilator_at_x_0=" + append.str() + "_with " + std::to_string(num_of_eigenVal) + " eigenvalues.txt";

		saveFileDialog.FileName = gcnew String(name.c_str());

		std::string name_base = "E_";
		std::string vec_base = "EigenVec_";

		if (saveFileDialog.ShowDialog() == System::Windows::Forms::DialogResult::OK) {
			System::IO::Stream^ strim = saveFileDialog.OpenFile();
			System::IO::StreamWriter^ outfile = gcnew System::IO::StreamWriter(strim);
			std::ostringstream energy_data, parameters;
			outfile->Write("-------------------------------------PARAMETERS: -------------------------------------\n");
			outfile->Write(L"Potential: Parabolic - Harmonic Oscilator\nV(x) = 1/2*m*omega^2*x^2\n\n"); 
			parameters << "Angular Frequency = " << std::setprecision(2) << omega << " [pHz]" << std::endl;
			parameters << "Quantum Well position = " << std::setprecision(2) << x_0 << " [nm]" << std::endl;
			parameters << "Quantum Well width = " << std::setprecision(2) << width << " [nm]" << std::endl;
			parameters << "x-axis sampling = = " << std::setprecision(4) << dx << " [nm]" << std::endl;
			outfile->Write(gcnew String((parameters.str()).c_str()));

			outfile->Write("\n-------------------------------------ENERGY DATA: -------------------------------------\
	\n\nno. state\tE_num [eV]\tE_analityc [eV]\n");
			for (int k = 0; k < num_of_eigenVal; k++) {
				energy_data << k << "\t\t" << std::setprecision(5) << eigenval_num[k] << "\t\t" << std::setprecision(5) << eigenVal_analityc[k] << std::endl;
			}
			outfile->Write(gcnew String((energy_data.str()).c_str()));
			if (this->checkBox_num->Checked) {
				outfile->Write("\n\n-------------------------------------Wavefunction Numerical: -------------------------------------\n\n");
				std::string	data = "x [nm]\t\t";
				std::ostringstream  wavefunction_data;
				for (int k = 0; k < num_of_eigenVal; k++) 
					data += std::to_string(k) + "-th eigenvector\t";
				wavefunction_data << data << "\n";
				for (int n = 0; n < N; n++) {
					double x = x_1 + n * dx;
					wavefunction_data << std::setprecision(4) << x << "\t\t";
					for (int k = 0; k < num_of_eigenVal; k++) {
						wavefunction_data << std::setprecision(10) << eigenVec_num[k][n] << "\t\t";
					}
					wavefunction_data << std::endl;
				}
				outfile->Write(gcnew String((wavefunction_data.str()).c_str()));
			}
			if (this->checkBox_analityc->Checked) {
				outfile->Write("\n\n-------------------------------------Wavefunction Analitycal: -------------------------------------\n\n");
				std::string	data = "x [nm]\t\t";
				std::ostringstream  wavefunction_data;
				for (int k = 0; k < num_of_eigenVal; k++)
					data += std::to_string(k) + "-th eigenvector\t";
				wavefunction_data << data << "\n";
				for (int n = 0; n < N; n++) {
					double x = x_1 + n * dx;
					wavefunction_data << std::setprecision(4) << x << "\t\t";
					for (int k = 0; k < num_of_eigenVal; k++) {
						wavefunction_data << std::setprecision(10) << eigenVec_analitycal[k][n] << "\t\t";
					}
					wavefunction_data << std::endl;
				}
				outfile->Write(gcnew String((wavefunction_data.str()).c_str()));
			}
			outfile->Close();
		}

		SetCurrentDirectory(Npath);
	}
};
}
